# 2dSYMPOL

**SYMmetry-based prediction of POLarity in 2D bilayers**

A Python tool for identifying polar and non-polar stacking configurations in 2D bilayer systems using crystallographic symmetry analysis.

## Features

- **Real-time search**: Find AA, AB, and BA stacking configurations for any 2D material
- **Symmetry-based analysis**: Uses layer group symmetry operations to classify stackings
- **Polarization direction**: Identifies x, y, z, and general polarization directions
- **High-symmetry prioritization**: Prioritizes physically meaningful stackings (e.g., 1/3, 2/3 for TMDCs)
- **CIF generation**: Creates bilayer structure files for visualization and analysis
- **c2db integration**: Works with the Computational 2D Materials Database

## Quick Start

```bash
# Search for WS2 stacking configurations
python3 run_2dsympol.py search --uid 1WS2-1

# Search by chemical formula
python3 run_2dsympol.py search --formula MoS2 --auto-select

# Find only z-polarized stackings with custom grid
python3 run_2dsympol.py search --uid 1WS2-1 --polar-direction z --grid 60

# List available materials
python3 run_2dsympol.py list --formula WS2
```

## Installation

1. Clone the repository
2. Install dependencies: `numpy`, `sqlite3`
3. Obtain the c2db database file (`c2db.db`) separately from c2db maintainers
4. Place `c2db.db` in the project root directory

## Mathematical Framework

Tests symmetry preservation using: `(E + R)τ = n` where R is a symmetry operation, τ is the stacking vector, and n is an integer. Stackings that break inversion/mirror symmetries are classified as polar.

For WS2, 2dSYMPOL identifies:
- **AA stacking**: τ = [0.000, 0.000] (non-polar)
- **AB stacking**: τ = [0.333, 0.333] (z-polar)
- **BA stacking**: τ = [0.667, 0.667] (z-polar, opposite to AB)

## Applications

- van der Waals heterostructure design
- Ferroelectric bilayer prediction
- Stacking-dependent property analysis
- 2D material database generation

## Documentation

See `MANUAL.md` for detailed usage instructions and methodology.

## Acknowledgments

This tool uses the Computational 2D Materials Database (c2db). Please obtain the database separately and cite:
```
Haastrup, S., Strange, M., Pandey, M. et al. 
The Computational 2D Materials Database: high-throughput modeling and discovery of atomically thin crystals. 
2D Mater. 5, 042002 (2018).
```

## Author

**Amal Kishore**

## Citation

If you use 2dSYMPOL in your research, please cite:
```
[Citation information to be added]
```

## License

[License information to be added]