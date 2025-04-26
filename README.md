# py-closewat

# The py-closewat project analyzes water molecules in X-ray diffraction PDB structures
# with resolutions between 1.01 Å and 1.06 Å. It includes tools to count waters,
# analyze distributions, and a Python port of the original C algorithm (closewat).

# Repository Structure
- `closewat.c`: Original C program for water molecule analysis.
- `pyclosewat.py`: Python reimplementation of `closewat.c`.
- `fetch_and_count_pdb_info.py`: Fetches PDB IDs via RCSB API and counts water molecules.
- `analyze_pdb_files.py`: Reads local PDB files in `pdb_files/`, extracts resolution, Rfree,
  and water count, outputs `pdb_water_analysis.csv`.
- `analyze_water_distribution.py`: Analyzes `pdb_water_analysis.csv` to compute statistics,
  generate bins, and produce plots (`water_distribution_histogram.png`,
  `resolution_vs_water.png`, `rfree_vs_water.png`).
- `pdb_files/`: Directory containing all X-ray diffraction PDB files of interest (1.01–1.06 Å).
- `requirements.txt`: Python dependencies for API fetching.
- CSV outputs: `pdb_water_analysis.csv`, `water_analysis_with_bins.csv`.
- Plots: `water_distribution_histogram.png`, `resolution_vs_water.png`, `rfree_vs_water.png`.
- `README.md`: This project overview.

# Usage
1. Ensure `pdb_files/` contains the desired PDB entries.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Extract PDB metadata and water counts:
   ```bash
   python analyze_pdb_files.py
   ```
4. Analyze distribution and generate plots:
   ```bash
   python analyze_water_distribution.py
   ```
5. (Optional) Run Python port on a single PDB:
   ```bash
   python pyclosewat.py [options] <input.pdb> <output.pdb>
   ```

# License
MIT