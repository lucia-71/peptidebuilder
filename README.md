# PeptideBuilder

A computational tool for designing and optimizing peptide structures through fragment-based molecular building and energy calculations using machine learning interatomic potentials.

## Overview

PeptideBuilder is a Python-based system that combines fragment-based molecular assembly with machine learning energy predictions to generate optimized peptide structures. The tool allows users to:

- Define amino acid fragments and ligand binding sites
- Randomly place and rotate fragments within a specified binding region
- Calculate interaction energies using the FAIRChem UMA-OMOL model
- Identify optimal fragment poses with spatial constraints
- Combine best poses to generate complete peptide sequences
- Estimate linker requirements between fragments

## Features

- **Fragment Library**: Pre-defined amino acid fragments including arginine, lysine, aspartic acid, glutamic acid, isoleucine, leucine, serine, tryptophan, tyrosine, and phenylalanine
- **Ligand Support**: Includes UV filter molecules (octinoxate, octocrylene, oxybenzone, avobenzone)
- **ML-Based Energy Calculation**: Uses FAIRChem's UMA-OMOL pre-trained model for accurate force field predictions
- **Pose Optimization**: BFGS geometry optimization for fragment placement
- **Spatial Analysis**: Calculates distances between fragment centers and suggests linker placement
- **3D Visualization**: Integration with py3Dmol for interactive molecular visualization
- **Flexible Configuration**: Supports customizable vdW distances, number of placement attempts, and charge/spin states

## Installation

### Dependencies

The tool requires the following Python packages:

```bash
pip install py3Dmol
pip install fairchem-core
pip install ase
pip install numpy
pip install pandas
pip install torch
pip install deepchem
```

### Quick Setup

```bash
git clone https://github.com/lucia-71/peptidebuilder.git
cd peptidebuilder
```

## Usage

### Basic Workflow

```python
from peptide_builder import *

# 1. Define fragment library
fragments = define_fragments()

# 2. Calculate spatial extent for each fragment
for frag in fragments:
    get_fragment_cutoff(frag, 0.5)

# 3. Initialize ML calculator
device = "cuda" if torch.cuda.is_available() else "cpu"
predictor = pretrained_mlip.get_predict_unit("uma-s-1", device=device)
calculator = FAIRChemCalculator(predictor, task_name="omol")

# 4. Define binding site dimensions
max_values, min_values = get_binding_site_dims(octinoxate_data, 0.0)

# 5. Generate fragment poses
new_molecules, ies = grow_fragments(
    [octinoxate_data],
    fragments,
    number_tries=100,
    calculator=calculator,
    max_values=max_values,
    min_values=min_values
)

# 6. Save results
save_xyz_files(new_molecules, fragments, octinoxate_data)

# 7. Select best non-overlapping poses
best_poses = get_best_seperate_poses(ies, fragments, octinoxate_data, new_molecules)

# 8. Combine and analyze
combined_poses, combined_atoms, all_centers, all_names = combine_best_poses(
    octinoxate_data, fragments, new_molecules, best_poses
)

# 9. Determine peptide sequence with linkers
distances = distance_between_centers(all_centers, all_names)
# or for sequential neighbor ordering:
distances = test_centers(all_centers, all_names)
```

### Key Functions

#### Fragment Management
- `define_fragments()` - Returns library of amino acid fragments
- `get_fragment_cutoff(frag, vdw_frac)` - Calculates spatial extent of fragment
- `show_fragment(frag)` - Displays fragment atoms and coordinates

#### Geometry
- `get_binding_site_dims(molecule, vdw_distance)` - Determines binding region boundaries
- `get_frag_coordinates(new_origin, which_frag)` - Generates rotated/translated fragment coordinates
- `calc_distance(ref, atom)` - Calculates Euclidean distance between points

#### Sampling & Optimization
- `grow_fragments(molecules, frags, tries, calculator, max_vals, min_vals)` - Samples fragment poses
- `calc_frag_energy(molecule, frag, symbols, calculator, ...)` - Computes interaction energy
- `distance_too_short(frag_coords1, frag_coords2, threshold)` - Checks for steric clashes

#### Pose Selection
- `get_best_poses(ies, frags)` - Ranks poses by energy
- `get_best_seperate_poses(ies, frags, bs, molecules)` - Selects non-overlapping poses
- `combine_best_poses(frags, bs, molecules, best_poses)` - Merges best poses

#### Analysis & Output
- `save_xyz_files(molecules, frags, molecule_dict)` - Saves structure files
- `view_frag_pose(frag_idx, pose_idx, frags, bs)` - Visualizes individual pose
- `view_combined_poses(bs, frags, best_poses)` - Visualizes combined structure
- `distance_between_centers(centers, names)` - Analyzes inter-fragment spacing
- `test_centers(centers, names)` - Determines sequential fragment ordering

## Data Format

### Fragment Dictionary Structure

```python
fragment = {
    "num_atoms": int,           # Total number of atoms
    "name": str,                # Fragment identifier
    "atoms": [str, ...],        # Element symbols
    "coords": [[x, y, z], ...], # 3D coordinates (Angstroms)
    "charge": int,              # Net charge
    "spin": int,                # Spin multiplicity
    "size": float               # Spatial extent (calculated)
}
```

## Output Files

- `frag_files/<ligand>_w_<fragment><index>.xyz` - Individual fragment poses in XYZ format
- `frag_files/combined.xyz` - Combined best poses

## Parameters

### Energy Calculation
- **Energy Cutoff**: 500.0 kcal/mol - Maximum acceptable energy for initial screening
- **Optimization Threshold**: fmax=0.05 eV/Angstrom - BFGS convergence criterion

### Sampling
- **vdW Fraction**: 0.5 - Multiplier for fragment size cutoff
- **Placement Attempts**: Configurable (typically 100-1000)
- **Steric Threshold**: 1.4 Angstroms - Minimum allowed interatomic distance

### Binding Site
- **vdW Buffer**: Configurable distance added to binding region boundaries

## Supported Molecules

### Amino Acids
Arginine, Lysine, Aspartic Acid, Glutamic Acid, Isoleucine, Leucine, Serine, Tryptophan, Tyrosine, Phenylalanine

### UV Filters (Example Binding Sites)
- Octinoxate (44 atoms)
- Octocrylene (56 atoms)
- Oxybenzone (29 atoms)
- Avobenzone (45 atoms)

## Computational Requirements

- **GPU**: Recommended (CUDA-compatible for faster calculations)
- **RAM**: Minimum 4GB (8GB+ recommended)
- **CPU**: Multi-core processor beneficial for parallel processing

## File Organization

```
peptidebuilder/
├── README.md                    # This file
├── peptide_builder.py           # Main module with core functions
├── frag_grow.ipynb             # Jupyter notebook with example workflow
├── frag_files/                 # Output directory for structure files
├── xyz_files/                  # Storage for XYZ coordinate files
├── md_videos/                  # Directory for MD trajectory outputs
└── log_files/                  # Directory for log files
```

## Example Notebook

The `frag_grow.ipynb` Jupyter notebook provides a complete workflow example including:
1. Environment setup and package installation
2. Fragment definition and cutoff calculation
3. ML model initialization
4. Fragment sampling and energy calculation
5. Best pose selection
6. Structure visualization and analysis
7. Peptide sequence determination

## Tips for Use

1. **Fragment Sampling**: Start with 100 placement attempts and increase if insufficient valid poses are found
2. **Geometry Optimization**: Disable optimization (`opt_flag=False`) for speed during initial screening
3. **Energy Thresholds**: Adjust the 500 kcal/mol cutoff based on your binding site characteristics
4. **Steric Clashes**: Reduce the 1.4 Angstrom threshold if poses are too close; increase for looser constraints
5. **Linker Estimation**: Use `distance_between_centers()` for shortest path; use `test_centers()` for sequential ordering

## Limitations

- Designed for peptide/small molecule binding site decoration
- ML predictions valid within FAIRChem training domain
- Fragment orientations limited to random rotations (no directed sampling)
- Interaction energies exclude long-range electrostatics beyond ML cutoff

## Future Enhancements

- Support for custom fragment libraries
- Directed sampling algorithms (e.g., genetic algorithms)
- Integration with protein structure databases
- Machine learning model fine-tuning capability
- Molecular dynamics post-processing

## Citation

If you use PeptideBuilder in your research, please cite:

```
PeptideBuilder: Fragment-based peptide design using machine learning
[Your citation information here]
```

## License

[Specify your license here - e.g., MIT, Apache 2.0, etc.]

## Contact

For questions or issues, please visit the [GitHub repository](https://github.com/lucia-71/peptidebuilder) or open an issue.

## Acknowledgments

- FAIRChem team for the UMA-OMOL ML interatomic potential
- ASE (Atomic Simulation Environment) for geometry optimization
- py3Dmol for molecular visualization

---

**Last Updated**: June 2026
