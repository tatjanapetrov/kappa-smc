# Statistical Model Checking for rule-based models in the Kappa language

This repository contains a research prototype for statistical model checking of stochastic models written in the [Kappa](https://kappalanguage.org/) language. The method is described in details here: https://doi.org/10.1007/978-3-032-30870-2_5

Given a Kappa model and a bounded linear temporal logic (BLTL) property, the framework:

1. parses the BLTL formula and extracts its atomic propositions;
2. instruments the Kappa model with the observables and interventions required to evaluate those propositions;
3. generates stochastic trajectories using KaSim through Kappy;
4. monitors the BLTL property over each trajectory;
5. estimates its satisfaction probability from independent simulation runs;
6. reports a Wilson confidence interval for the estimated probability.

The repository includes examples from systems biology and collective decision-making.

## Repository structure

```text
kappa-smc/
├── pyproject.toml
├── README.md
├── src/
│   └── kappa_smc/
│       ├── __init__.py
│       ├── bltl.py                 # BLTL formula representation, parsing, and offline monitoring
│       ├── instrumentation.py      # automatic instrumentation of Kappa models
│       ├── simulation.py           # execution of KaSim simulations and processing of trajectory data
│       └── estimation.py           # statistical model checking and result structures
├── models/
├── experiments/
├── examples/
└── notebooks/
    └── demonstration.ipynb         # executable demonstration of the complete workflow
```

## Requirements

The implementation requires:

- Python 3.10 or later;
- NumPy;
- Matplotlib;
- Kappy;
- a local Kappa installation containing the KaSim executables.

Kappa must be installed separately. The demonstration notebook requires the path to the local Kappa executable directory.

## Installation

Clone the repository and enter its root directory:

```bash
git clone https://github.com/tatjanapetrov/kappa-smc
cd kappa-smc
```

Creating a dedicated virtual environment is recommended.

### Windows

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

### Linux and macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

The editable installation makes the `kappa_smc` package available while keeping the source files directly editable.

Verify the installation with:

```bash
python -c "import kappa_smc; print(kappa_smc.__file__)"
```

### Kappy on Windows

Kappy may not be installable directly from PyPI on native Windows. On Windows, a working Kappy installation must already be available in the selected Python environment.

Verify it with:

```powershell
python -c "import kappy; print(kappy.__file__)"
```

The notebook kernel and the terminal used for installation must use the same Python environment.

## Supported BLTL syntax

The current parser supports the following constructors:

```text
AP("variable>threshold")
Not(phi)
And(phi, psi)
Or(phi, psi)
F_leq(bound, phi)
G_leq(bound, phi)
Until_leq(bound, phi, psi)
```

Examples:

```python
'F_leq(50, AP("I<2"))'
```

```python
(
    "And("
    'G_leq(50, AP("I<80")),'
    'F_leq(50, AP("I<2"))'
    ")"
)
```

Atomic propositions currently consist of:

- the name of a Kappa algebraic variable;
- a strict comparison operator, `<` or `>`;
- an integer threshold.

## Cite us

> Salazar, A., Feret, J., & Petrov, T. (2027). Statistical Model Checking for Rule-Based Models in the Kappa Language. In P. T. Monteiro (Ed.), Computational Methods in Systems Biology: CMSB 2026 (Lecture Notes in Computer Science, Vol. 16678, pp. 73–103). Springer, Cham. https://doi.org/10.1007/978-3-032-30870-2_5

```bibtex
@inproceedings{SalazarFeretPetrov2027,
  author    = {Salazar, Albin and Feret, J{\'e}r{\^o}me and Petrov, Tatjana},
  editor    = {Monteiro, Pedro T.},
  title     = {Statistical Model Checking for Rule-Based Models in the Kappa Language},
  booktitle = {Computational Methods in Systems Biology},
  series    = {Lecture Notes in Computer Science},
  volume    = {16678},
  pages     = {73--103},
  year      = {2027},
  publisher = {Springer Nature Switzerland},
  address   = {Cham},
  doi       = {10.1007/978-3-032-30870-2_5},
  isbn      = {978-3-032-30870-2}
}
```

## Licence

See the `LICENSE` file for the terms under which this software is distributed.