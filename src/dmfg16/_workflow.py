import sys
import os
import subprocess
import tempfile
from pathlib import Path
from ._io import parse_gaussian_qst_input,sanitize_route,convert_gaussian_qst_in
from dmf import DirectMaxFlux, interpolate_fbenm
from ase.calculators.gaussian import Gaussian
from ase.io import read, write
import numpy as np


def read_input(inp_qst):

    with open(inp_qst) as fd:
        cgs = parse_gaussian_qst_input(fd)

    params = cgs[0].parameters.copy()

    is_qst = False
    for k,v in params.items():
        if ((isinstance(k,str) and "qst" in k)
            or (isinstance(v,str) and "qst" in v)):
            is_qst = True
            break

    if not is_qst:
        return is_qst, None, None

    if len(cgs)<2:
        raise ValueError(
            f"PyDMF: QST input requires multiple structures, but got {len(cgs)}"
        )

    ref_images = [cgs[0].get_atoms()]
    for cg in cgs[2:]:
        ref_images.append(cg.get_atoms())
    ref_images.append(cgs[1].get_atoms())

    return is_qst, ref_images, params


def get_coefs(ref_images,p_dir_dmf):

    images_ext = [ref_images[0]]
    nrefs = len(ref_images)

    if nrefs>2:
        sequential = False
    else:
        sequential = True

    for i in range(nrefs-1):

        mxflx_fbenm = interpolate_fbenm(
            [images_ext[-1],ref_images[i+1]],
            sequential = sequential,
            output_file = str(p_dir_dmf / f"ipopt_fbenm{i}.log"),
            ipopt_options = {
                "print_level": 0,
                "file_print_level": 5,
            },

        )

        if nrefs==2:
            return mxflx_fbenm.coefs

        images_ext.extend(mxflx_fbenm.images[1:])

    mxflx = DirectMaxFlux(images_ext)
    return mxflx.coefs


def get_ts_guess(ref_images,params,args):

    nmove = args.npoints
    p_dir_dmf = Path(args.dir)
    exe = args.exe
    parallel = args.parallel
    tol = args.tol
    restart = args.restart
    update_teval = not args.equidistant

    p_coefs = p_dir_dmf / "coefs.npy"
    p_image00 = p_dir_dmf / "image00.log"

    if restart and p_coefs.exists():
        coefs = np.load(p_coefs)
    elif restart and p_image00.exists():
        rst_images = []
        iimg = 0
        p_next = p_dir_dmf / f"image{iimg:02}.log"
        while p_next.exists():
            rst_images.append(read(p_next,format="gaussian-out"))
            iimg += 1
            p_next = p_dir_dmf / f"image{iimg:02}.log"
        coefs = get_coefs(rst_images,p_dir_dmf)
    else:
        coefs = get_coefs(ref_images,p_dir_dmf)

    params4calc = sanitize_route(params)
    params4calc.pop("chk",None)
    params4calc.pop("oldchk",None)

    if parallel:
        parallel, cpumems, params4calc = check_parallel(params4calc,nmove)

    mxflx = DirectMaxFlux(
        ref_images,
        coefs=coefs,
        nmove=nmove,
        update_teval=update_teval,
    )

    ipopt_options = {
        "print_level": 0,
        "output_file": str(p_dir_dmf / "ipopt_dmf.log"),
        "file_print_level": 5,
    }

    mxflx.add_ipopt_options(ipopt_options)
    command = exe+" < PREFIX.com > PREFIX.log"

    for i,image in enumerate(mxflx.images):
        label = str(p_dir_dmf / f"image{i:02}")
        image.calc = Gaussian(
            command=command,
            label=label,
            chk=f"image{i:02}.chk",
            ioplist=['2/11=1','2/12=1'],
            **params4calc
        )

    if parallel:
        mxflx.parallel = "thread"
        for image,cm in zip(mxflx.images,cpumems):
            image.calc.set(**cm)

    mxflx.get_forces()
    for image in mxflx.images:
        image.calc.set(**{"guess":"read"})

    mxflx.solve(tol=tol)

    np.save(str(p_dir_dmf / "coefs"), mxflx.coefs)
    np.save(str(p_dir_dmf / "tevals"), mxflx.coefs)
    write(str(p_dir_dmf / "dmf.traj"), mxflx.images)

    return mxflx.history.images_tmax[-1]


def check_parallel(params4calc,nmove):

    params4calc_ini = params4calc.copy()

    is_gpu = False
    is_linda = False

    if os.environ.get("GAUSS_GDEF") is not None:
        is_gpu = True

    for k in params4calc:
        if "gpu" in k.lower():
            is_gpu = True
        if "linda" in k.lower():
            is_linda = True

    if is_gpu:
        print(
           "PyDMF: DMF-parallel disabled due to Gaussian GPU calculation."
        )
        return False, None, params4calc_ini

    if is_linda:
        print(
            "PyDMF: DMF-parallel disabled due to Gaussian lindaworkers."
        )
        return False, None, params4calc_ini

    cpu = params4calc.pop("cpu",None)
    if cpu is None:
        cpu = os.environ.get("GAUSS_CDEF")
    if cpu is None:
        nprocshared = params4calc.pop("nprocshared",None)
    if cpu is None and nprocshared is None:
        print(
            "PyDMF: DMF-parallel disabled. Either cpu or nprocshared must be set in Link0 section or environment variables."
        )
        return False, None, params4calc_ini

    if cpu is not None:
        cpu_split = cpu.split(",")
        cores = []
        for c in cpu_split:
            c_split = list(map(int,c.split("-")))
            if len(c_split)==1:
                cores.append(c_split[0])
            else:
                cores.extend(list(map(str,range(c_split[0],c_split[1]+1))))
        cores = list(map(str,cores))
        nproc = len(cores)
    else:
        nproc = int(nprocshared)

    if nproc < nmove:
        print(
            "PyDMF: DMF-parallel disabled. Number of CPU cores lower than DMF movable points."
        )
        return False, None, params4calc_ini

    cpumems = []
    nproc_per_img = nproc // nmove
    if cpu is not None:
        cpumems.append({"cpu":",".join(cores[:nproc//2])})
        for i in range(nmove):
            cpumems.append(
                {"cpu":",".join(cores[i*nproc_per_img:(i+1)*nproc_per_img])}
            )
        cpumems.append({"cpu":",".join(cores[nproc//2:])})
    else:
        cpumems.append({"nprocshared":str(nproc//2)})
        for i in range(nmove):
            cpumems.append({"nprocshared":str(nproc_per_img)})
        cpumems.append({"nprocshared":str(nproc//2)})

    mem =  params4calc.pop("mem",None)
    if mem is not None:
        mem = mem.strip(" ")
        mem_val = int(mem[:-2])
        mem_unit = mem[-2:]

        if mem_val<100:
            kmgt = "kmgt"
            for i in range(3):
                if mem_unit[0]==kmgt[i+1]:
                    mem_val *= 1000
                    mem_unit = kmgt[i]+mem_unit[1]
                    break

        for i, cm in enumerate(cpumems):
            if i==0 or i==len(cpumems)-1:
                cm["mem"] = f"{mem_val//2}{mem_unit}"
            else:
                cm["mem"] = f"{mem_val//nmove}{mem_unit}"

    return True, cpumems, params4calc


def main(args):

    nmove = args.npoints
    dir_dmf = args.dir
    exe = args.exe
    parallel = args.parallel
    only_ts_guess = args.only_ts_guess

    inp = sys.stdin.read()
    if not inp.strip():
        raise RuntimeError("PyDMF: No input provided on stdin")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        com = tmpdir / "input.com"
        log = tmpdir / "output.log"

        com.write_text(inp)

        is_qst, ref_images, params = read_input(com)

        if is_qst:
            p_dir_dmf = Path(dir_dmf)
            p_dir_dmf.mkdir(exist_ok=True)

            atoms_guess = get_ts_guess(ref_images,params,args)
            new_inp = convert_gaussian_qst_in(com,atoms_guess,params)
            new_com = p_dir_dmf / "converted.com"
            new_com.write_text(new_inp)

            print("PyDMF: QST option found. Converted input is submitted.")
            com2gau = new_com
        else:
            print("PyDMF: QST option not found. Original input is submitted.")
            com2gau = com

        if only_ts_guess:
            print("PyDMF: Gaussian excusion skipped.")
            return

        cmd = [exe]

        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=None,
            bufsize=1,
            text=True,
        )

        with open(com2gau, "r") as fin:
            proc.stdin.write(fin.read())
        proc.stdin.close()

        with open(log, "w") as fout:
            for line in proc.stdout:
                fout.write(line)
                fout.flush()
                sys.stdout.write(line)
                sys.stdout.flush()

        proc.wait()

        sys.exit(proc.returncode)


