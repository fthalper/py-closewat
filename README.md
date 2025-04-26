# py-closewat

[![Python Version](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Python toolkit to analyze water molecules in X-ray diffraction PDB structures with resolutions between 1.01 Å and 1.06 Å. It provides utilities to count waters, analyze distributions, and includes a Python reimplementation of the original C `closewat` algorithm.

## Table of Contents

- [Features](#features)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [Outputs](#outputs)
- [Development](#development)
- [License](#license)
- [Citation](#citation)

## Features

- Count water molecules in PDB files
- Extract resolution and R-free values
- Bin and analyze water distributions across high-resolution structures
- Generate plots: histograms, resolution vs water count, R-free vs water count
- Python and C implementations of the `closewat` algorithm

## Getting Started

### Prerequisites

- Python 3.7+
- `pip`

### Installation

Install the required Python packages:

```bash
pip install -r requirements.txt
```

## Usage

### Analyze Multiple PDB Structures

Place your PDB files (1.01–1.06 Å) in the `pdb_files/` directory, then run:

```bash
python analyze_pdb_files.py
python analyze_water_distribution.py
```

This will produce:
- `pdb_water_analysis.csv`
- `water_analysis_with_bins.csv`
- `water_distribution_histogram.png`
- `resolution_vs_water.png`
- `rfree_vs_water.png`

### Analyze a Single PDB File

Using the Python port:

```bash
python pyclosewat.py [options] <input.pdb> <output.pdb>
```

Example:

```bash
python pyclosewat.py -d 3.5 1abc.pdb 1abc_annotated.pdb
```

Using the C program:

```bash
gcc closewat.c -o closewat
./closewat -d 3.5 1abc.pdb 1abc_annotated.pdb
```

## Repository Structure

```
py-closewat/
├── closewat.c
├── pyclosewat.py
├── fetch_and_count_pdb_info.py
├── analyze_pdb_files.py
├── analyze_water_distribution.py
├── pdb_files/
├── requirements.txt
├── README.md
└── LICENSE
```

## Outputs

- **pdb_water_analysis.csv**: Raw PDB metadata and water counts
- **water_analysis_with_bins.csv**: Binned water count data
- **Plots**: PNG files for distribution and correlation analyses

## Development

Contributions welcome! Please open issues or submit pull requests.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use this toolbox in your research, please cite:

> Author Name. _py-closewat_: analyzing water molecules in PDB structures. GitHub repository. Version X.X, Year.
