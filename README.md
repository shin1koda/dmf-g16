# dmf-g16: A Gaussian wrapper for PyDMF double-ended transition-state searches

## Requirements

- [ASE](https://ase-lib.org/)
- [cyipopt](https://cyipopt.readthedocs.io/en/stable/)
- [PyDMF](https://github.com/shin1koda/dmf)


## Installation

We generally recommend installing this package via **conda**, as `cyipopt` is most reliably installed through conda.

```bash
conda create -n dmfg16 python=3.10
conda activate dmfg16
conda install -c conda-forge ase cyipopt
pip install dmfg16
```


## Usage

Just replace excutable from `g16` to `dmf-g16` as follows.

```bash
#g16 < input.com > log
dmf-g16 < input.com > log
```


## Citation

 1. S.-i. Koda and  S. Saito, Locating Transition States by Variational Reaction Path Optimization with an Energy-Derivative-Free Objective Function, JCTC, 20, 2798–2811 (2024). [doi: 10.1021/acs.jctc.3c01246](https://doi.org/10.1021/acs.jctc.3c01246)
 1. S.-i. Koda and  S. Saito, Flat-bottom Elastic Network Model for Generating Improved Plausible Reaction Paths, JCTC, 20, 7176−7187 (2024). [doi: 10.1021/acs.jctc.4c00792](https://doi.org/10.1021/acs.jctc.4c00792)
 1. S.-i. Koda and  S. Saito, Correlated Flat-bottom Elastic Network Model for Improved Bond Rearrangement in Reaction Paths, JCTC, 21, 3513−3522 (2025). [doi: 10.1021/acs.jctc.4c01549](https://doi.org/10.1021/acs.jctc.4c01549)


## Community guidelines

### Contributing

Contributions to this project are welcome. If you would like to contribute new features, improvements, or documentation, please open a pull request on GitHub.  
Before submitting a PR, we recommend opening a short issue to discuss the proposed change.

### Reporting issues

If you encounter a problem, unexpected behavior, or a potential bug, please report it through the GitHub issue tracker:

https://github.com/shin1koda/dmf-g16/issues

When reporting an issue, please include:
- A clear description of the problem  
- Steps to reproduce the issue  
- Your environment (Python version, ASE version, cyipopt version, etc.)  
- Any relevant error messages or logs

### Seeking support

If you have questions about the usage of the package, or need help integrating it into your workflow, feel free to open an issue labeled “question” on GitHub.  
We will do our best to provide guidance based on availability.


## License

This software is licensed under the GNU Lesser General Public License v2.1 or later.
This software includes modified code derived from the Atomic Simulation Environment (ASE).
