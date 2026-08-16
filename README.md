# dmf-g16: A Gaussian wrapper for PyDMF double-ended transition-state searches

## Overview

dmf-g16 is a Gaussian wrapper for Direct MaxFlux (DMF)-based double-ended transition-state (TS) searches. It allows Gaussian users to perform DMF-based reaction-path optimization through PyDMF while keeping their existing Gaussian workflows almost unchanged.

Users can run dmf-g16 with native Gaussian QST2/QST3 input files by simply replacing the Gaussian executable, such as `g16`, with `dmf-g16`. For QST inputs, dmf-g16 performs DMF-based path optimization using Gaussian for energy and gradient evaluations, then runs a Gaussian TS optimization from the highest-energy point on the optimized path.

## Platform support

dmf-g16 supports Linux and Windows environments.

We gratefully acknowledge Dr. Hideya Tanaka (@tanaka-hideya) for contributing Windows support to dmf-g16.

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

If you use dmf-g16 in your research, please cite Refs. 1 and 2.

 1. S.-i. Koda and  S. Saito, dmf-g16: A Gaussian Wrapper for Reliable Double-Ended Transition-State Searches With Native Input Formats, JCC, 47, e70378 (2026). [doi: 10.1002/jcc.70378](https://doi.org/10.1002/jcc.70378)
 2. S.-i. Koda and  S. Saito, PyDMF: A Python package for variational double-ended reaction-path and transition-state optimization, J. Open Source Softw. 11, 10379 (2026). [doi: 10.21105/joss.10379](https://doi.org/10.21105/joss.10379)


The underlying methods used by dmf-g16 are described in Refs. 3–5. We would greatly appreciate it if you could also cite these papers.

 3. S.-i. Koda and  S. Saito, Locating Transition States by Variational Reaction Path Optimization with an Energy-Derivative-Free Objective Function, JCTC, 20, 2798–2811 (2024). [doi: 10.1021/acs.jctc.3c01246](https://doi.org/10.1021/acs.jctc.3c01246)
 4. S.-i. Koda and  S. Saito, Flat-bottom Elastic Network Model for Generating Improved Plausible Reaction Paths, JCTC, 20, 7176−7187 (2024). [doi: 10.1021/acs.jctc.4c00792](https://doi.org/10.1021/acs.jctc.4c00792)
 5. S.-i. Koda and  S. Saito, Correlated Flat-bottom Elastic Network Model for Improved Bond Rearrangement in Reaction Paths, JCTC, 21, 3513−3522 (2025). [doi: 10.1021/acs.jctc.4c01549](https://doi.org/10.1021/acs.jctc.4c01549)

For example:

> Double-ended transition-state searches were performed using dmf-g16 [1], a Gaussian [appropriate reference] wrapper for PyDMF [2], employing the Direct MaxFlux method [3] with two potential-energy models: the flat-bottom elastic network model [4,5] implemented in PyDMF for initial-path generation and DFT evaluated with Gaussian for obtaining approximate transition-state structures. These structures were further refined using Gaussian. All DFT calculations were performed at ... (description of the calculation level).

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
