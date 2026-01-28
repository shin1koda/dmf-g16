# This file includes code derived from the Atomic Simulation Environment (ASE)
# Copyright (C) ASE developers

import ase.io.gaussian as iog
import re
from pathlib import Path
from copy import deepcopy


def convert_gaussian_qst_in(inp,atoms,parameters):

    with open(inp) as f:
        out = _get_gaussian_converted_qst_in(
            f,
            atoms,
            **parameters
        )

    return "".join(out)

def _get_gaussian_converted_qst_in(fd,atoms,
    method=None, basis=None, fitting_basis=None,
    output_type='P', basisfile=None, basis_set=None,
    xc=None, charge=None, mult=None, extra=None,
    ioplist=None, addsec=None, spinlist=None,
    zefflist=None, qmomlist=None, nmagmlist=None,
    znuclist=None, radnuclearlist=None,
    **parameters):
    # These variables indicate which section of the
    # input file is currently being read:
    route_section = False
    atoms_section = False
    atoms_saved = False
    after_link1 = False
    # Here we will store the sections of the file in a dictionary,
    # as lists of strings
    gaussian_sections = {'link0': [], 'route': [],
                         'charge_mult': [], 'mol_spec': [], 'extra': [],
                         'link1': [],}

    idx_chgmult = -1

    for line in fd:
        if 'link1' in line.lower():
            after_link1 = True
        #line = line.strip(' ')
        link0_match = iog._re_link0.match(line)
        output_type_match = iog._re_output_type.match(line)
        chgmult_match = iog._re_chgmult.match(line)
        if after_link1:
            gaussian_sections["link1"].append(line)
        elif link0_match:
            gaussian_sections['link0'].append(line)
        elif line.strip(' ') == '\n' and (route_section or atoms_section):
            route_section = False
            atoms_section = False
        elif output_type_match or route_section:
            route_section = True
            gaussian_sections['route'].append(
                re.sub(
                    r'\bqst[23]\b',
                    'ts',
                    line,
                    flags=re.IGNORECASE,
                )
            )
        elif chgmult_match:
            idx_chgmult += 1
            gaussian_sections['charge_mult'] = line
            gaussian_sections[f'mol_spec{idx_chgmult}'] = []
            gaussian_sections['extra'] = []
            atoms_section = True
        elif atoms_section:
            gaussian_sections[f'mol_spec{idx_chgmult}'].append(line)
            atoms_saved = True
        elif atoms_saved:
            gaussian_sections['extra'].append(line)

    for k in ['mol_spec']:
        if f'{k}0' in gaussian_sections:
            gaussian_sections[k] = gaussian_sections[f'{k}0']

    # header, charge, and mult
    new_coords = ['\n', 'Gaussian input prepared by ASE and PyDMF\n', '\n',
                  f'{charge:.0f} {mult:.0f}\n']

    # make dict of nuclear properties:
    nuclear_props = {'spin': spinlist, 'zeff': zefflist, 'qmom': qmomlist,
                     'nmagm': nmagmlist, 'znuc': znuclist,
                     'radnuclear': radnuclearlist}
    nuclear_props = {k: v for k, v in nuclear_props.items() if v is not None}

    # atomic positions and nuclear properties:
    molecule_spec = iog._get_molecule_spec(atoms, nuclear_props)
    for line in molecule_spec:
        new_coords.append(line+"\n")

    out = []
    out.extend(gaussian_sections["link0"])
    out.extend(gaussian_sections["route"])
    out.extend(new_coords)
    out.extend(gaussian_sections["extra"])
    out.extend(gaussian_sections["link1"])

    return out

def _get_gaussian_qst_in_sections(fd):
    ''' Reads a gaussian input file and returns
    a dictionary of the sections of the file - each dictionary
    value is a list of strings - each string is a line in that
    section. '''
    # These variables indicate which section of the
    # input file is currently being read:
    route_section = False
    atoms_section = False
    atoms_saved = False
    # Here we will store the sections of the file in a dictionary,
    # as lists of strings
    gaussian_sections = {'link0': [], 'route': [],
                         'charge_mult': [], 'mol_spec': [], 'extra': []}

    idx_chgmult = -1

    for line in fd:
        if 'link1' in line.lower():
            break
        line = line.strip(' ')
        link0_match = iog._re_link0.match(line)
        output_type_match = iog._re_output_type.match(line)
        chgmult_match = iog._re_chgmult.match(line)
        if link0_match:
            gaussian_sections['link0'].append(line)
        # The first blank line appears at the end of the route section
        # and a blank line appears at the end of the atoms section:
        elif line == '\n' and (route_section or atoms_section):
            route_section = False
            atoms_section = False
        elif output_type_match or route_section:
            route_section = True
            gaussian_sections['route'].append(line)
        elif chgmult_match:
            idx_chgmult += 1
            gaussian_sections[f'charge_mult'] = line
            gaussian_sections[f'mol_spec{idx_chgmult}'] = []
            gaussian_sections[f'extra'] = []
            # After the charge and multiplicty have been set, the
            # molecule specification section of the input file begins:
            atoms_section = True
        elif atoms_section:
            gaussian_sections[f'mol_spec{idx_chgmult}'].append(line)
            atoms_saved = True
        elif atoms_saved:
            # Next we read the other sections below the molecule spec.
            # This may include the ReadIso section and the basis set
            # definition or filename
            gaussian_sections[f'extra'].append(line)

    for k in ['mol_spec']:
        if f'{k}0' in gaussian_sections:
            gaussian_sections[k] = gaussian_sections[f'{k}0']

    return gaussian_sections


def parse_gaussian_qst_input(fd):

    file_sections = _get_gaussian_qst_in_sections(fd)

    nmols = -1
    for k in file_sections:
        if "mol_spec" in k:
            nmols += 1

    gcs = []
    for idx in range(nmols):
        for k in ['mol_spec']:
            if f'{k}{idx}' in file_sections:
                file_sections[k] = file_sections[f'{k}{idx}']

        parameters = {}
        # Update the parameters dictionary with the keywords and values
        # from each section of the input file:
        parameters.update(iog._get_all_link0_params(file_sections['link0']))
        parameters.update(iog._get_all_route_params(file_sections['route']))
        parameters.update(iog._get_charge_mult(file_sections['charge_mult']))
        atoms, nuclear_props = iog._get_atoms_from_molspec(
            file_sections['mol_spec'])
        parameters.update(nuclear_props)
        parameters = _get_extra_section_params(
            file_sections['extra'], parameters, atoms)

        if parameters.get("basis_set","") != '':
            parameters.pop("basis",None)

        iog._validate_params(parameters)

        gcs.append(iog.GaussianConfiguration(atoms, parameters))

    return gcs


def read_gaussian_qst_in(fd, attach_calculator=False):

    gcs = parse_gaussian_qst_input(fd)

    list_atoms = []
    for gc in gcs:
        atoms = gc.get_atoms()

        if attach_calculator:
            atoms.calc = gc.get_calculator()

        list_atoms.append(atoms)

    return list_atoms

def sanitize_route(params):

    FORCE_INCOMPATIBLE_KEYS = {
        "sp", "opt", "freq", "irc", "ircmax", "scan",
        "polar", "admp", "bomd", "eet", "force",
        "stable", "volume", "geom", "modredundant",
        "ts", "qst2", "qst3", "calcfc", "calcall",
        "addredundant", "maxcycles", "maxstep", "optcyc",
        "restart", "readfc", "rcfc", "readcartesianfc",
        "recalcfc", "noraman", "fccards", "eigentest",
        "noeigentest", "gediis", "rfo", "ef", "eigenvaluefollow",
        "gic", "path", "raman", "nraman", "nnraman",
        "forward", "reverse", "maxpoints", "stepsize",
    }

    out = {
        k:v
        for k,v in params.items()
        if k.lower() not in FORCE_INCOMPATIBLE_KEYS
    }

    return out

def _get_extra_section_params(extra_section, parameters, atoms):
    ''' Takes a list of strings: extra_section, which contains
    the 'extra' lines in a gaussian input file. Also takes the parameters
    that have been read so far, and the atoms that have been read from the
    file.
    Returns an updated version of the parameters dict, containing details from
    this extra section. This may include the basis set definition or filename,
    and/or the readiso section.'''

    new_parameters = deepcopy(parameters)
    # Will store the basis set definition (if set):
    basis_set = ""
    lines_basis_set = []
    # This will indicate whether we have a readiso section:
    readiso = iog._get_readiso_param(new_parameters)[0]
    # This will indicate which line of the readiso section we are reading:
    count_iso = 0
    readiso_masses = []

    for line in extra_section:
        if line.split():
            # check that the line isn't just a comment
            if line.split()[0] == '!':
                continue
            line = line.strip().split('!')[0]

        #if len(line) > 0 and line[0] == '@':
        #    # If the name of a basis file is specified, this line
        #    # begins with a '@'
        #    new_parameters['basisfile'] = line
        #elif readiso and count_iso < len(atoms.symbols) + 1:
        if readiso and count_iso < len(atoms.symbols) + 1:
            # The ReadIso section has 1 more line than the number of
            # symbols
            if count_iso == 0 and line != '\n':
                # The first line in the readiso section contains the
                # temperature, pressure, scale. Here we save this:
                readiso_info = iog._get_readiso_info(line, new_parameters)
                if readiso_info is not None:
                    new_parameters.update(readiso_info)
                # If the atom masses were set in the nuclear properties
                # section, they will be overwritten by the ReadIso
                # section
                readiso_masses = []
                count_iso += 1
            elif count_iso > 0:
                # The remaining lines in the ReadIso section are
                # the masses of the atoms
                try:
                    readiso_masses.append(float(line))
                except ValueError:
                    readiso_masses.append(None)
                count_iso += 1
        else:
            # If the rest of the file is not the ReadIso section,
            # then it must be the definition of the basis set.
            if new_parameters.get('basis', '') == 'gen' \
                    or 'gen' in new_parameters:
                lines_basis_set.append(line.strip())

    if readiso:
        new_parameters['isolist'] = readiso_masses
        new_parameters = iog._update_readiso_params(new_parameters, atoms.symbols)

    basis_set = "\n".join(lines_basis_set).strip()

    # Saves the basis set definition to the parameters array if
    # it has been set:
    if basis_set != '':
        new_parameters['basis_set'] = basis_set
        new_parameters['basis'] = 'gen'
        new_parameters.pop('gen', None)

    return new_parameters

