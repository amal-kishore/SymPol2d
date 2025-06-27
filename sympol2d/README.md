# SYMPOL2D: SYMmetry-based prediction of POLarity in 2D bilayers

## Complete CLI Command Reference

### Main Command
```bash
sympol2d [--database DATABASE] {search,list} [options]
```

### Global Options
- `--database, -d DATABASE`: Path to C2DB database file (default: c2db.db)
- `--help, -h`: Show help message

### Search Command
```bash
sympol2d search [options]
```

#### Search Options:
- `--uid, -u UID`: Material UID (e.g., MoS2-165798ab5e18)
- `--formula, -f FORMULA`: Chemical formula (e.g., MoS2)
- `--grid, -g GRID`: Grid size for scanning (default: 50)
- `--output, -o OUTPUT`: Output path (directory for CIF, file for JSON)
- `--output-format, -of {cif,json}`: Output format (default: cif)
- `--auto-select`: Auto-select first match when multiple materials found
- `--polar-direction, -pd {x,y,z,xy,general,all}`: Filter by polarization direction (default: all)
- `--interlayer-distance, -d DISTANCE`: Interlayer distance in Angstroms (default: auto-estimate)

### List Command
```bash
sympol2d list [options]
```

#### List Options:
- `--formula, -f FORMULA`: Filter by chemical formula
- `--layer-group, -lg GROUP`: Filter by layer group
- `--layer-groups`: List all layer groups with material counts
- `--limit, -n LIMIT`: Maximum number of results (default: 20)

## Quick Start Examples

### Basic Usage

```bash
# Search for polar stackings by material UID (saves CIF files by default)
sympol2d search --uid MoS2-165798ab5e18

# Search by chemical formula (auto-select first match)
sympol2d search --formula MoS2 --auto-select

# Specify output directory for CIF files
sympol2d search --formula MoS2 --auto-select --output MoS2_bilayers

# Save results as JSON instead of CIF
sympol2d search --formula MoS2 --auto-select --output-format json --output results.json

# List materials in database
sympol2d list --formula WS2

# List all layer groups
sympol2d list --layer-groups
```

### Advanced Usage

```bash
# High-resolution search with specific polarization directions
sympol2d search --formula MoS2 --auto-select --grid 100 --polar-direction x
sympol2d search --uid 1WS2-1 --grid 80 --polar-direction z
sympol2d search --formula WSe2 --auto-select --polar-direction y
sympol2d search --formula MoTe2 --auto-select --polar-direction xy

# Using python script directly (if not installed)
python3 run_sympol2d.py search --uid 1WS2-1 --grid 80 --polar-direction z
python3 run_sympol2d.py search --formula MoS2 --auto-select --polar-direction x

# Custom interlayer distance and save CIF files
sympol2d search --formula MoS2 --auto-select --interlayer-distance 3.5 --output MoS2_d3.5
sympol2d search --uid 1WS2-1 --interlayer-distance 3.2 --output WS2_d3.2

# Save as JSON for data analysis
sympol2d search --formula MoS2 --auto-select --output-format json --output MoS2_results.json
sympol2d search --uid 1WS2-1 --grid 150 --output-format json --output WS2_fine_grid.json

# Fine grid scanning with CIF output
sympol2d search --formula MoS2 --auto-select --grid 100 --output MoS2_fine
sympol2d search --uid 1WS2-1 --grid 150 --polar-direction all --output WS2_all_polar

# Search with custom database
sympol2d --database /path/to/custom.db search --formula WS2 --auto-select
sympol2d --database ~/databases/c2db_latest.db search --uid 1WS2-1
```

## Handling Multiple Polar Pairs

When SYMPOL2D finds multiple polar pairs, it prioritizes them by symmetry:

1. **High-symmetry pairs** (e.g., τ = [1/3, 1/3] and [2/3, 2/3]) are listed first
2. **Lower-symmetry pairs** follow in order of decreasing symmetry

### What to do with multiple pairs:

1. **For initial screening**: Use the first (highest-symmetry) pair, which is usually the most stable configuration

2. **For comprehensive analysis**: 
   - Save all results to JSON: `--output all_pairs.json`
   - Each pair represents a different ferroelectric configuration
   - Different pairs may have different polarization strengths and switching barriers

3. **For specific applications**:
   - **Memory devices**: Choose pairs with large polarization difference
   - **Sensors**: Select pairs with specific polar directions (use `--polar-direction`)
   - **Photovoltaics**: Consider all pairs as they may have different band alignments

### Example: Working with multiple pairs

```bash
# Find all polar pairs for different materials
sympol2d search --formula MoS2 --auto-select --output MoS2_all_pairs.json
sympol2d search --uid 1WS2-1 --output WS2_all_pairs.json

# Find only specific polarization directions
sympol2d search --formula MoS2 --auto-select --polar-direction x
sympol2d search --uid 1WS2-1 --polar-direction z --grid 100
sympol2d search --formula WSe2 --auto-select --polar-direction y
sympol2d search --formula MoTe2 --auto-select --polar-direction xy

# Systematic study with different parameters
sympol2d search --formula MoS2 --auto-select --grid 50 --polar-direction all --output MoS2_g50.json
sympol2d search --formula MoS2 --auto-select --grid 100 --polar-direction all --output MoS2_g100.json
sympol2d search --formula MoS2 --auto-select --grid 150 --polar-direction all --output MoS2_g150.json

# The output will show:
# - Pair 1: Highest symmetry (e.g., standard AB/BA stacking)
# - Pair 2, 3, ...: Alternative polar configurations
```

### Understanding the output:

```
X-POLARIZED PAIRS (12 found):

  Pair 1:                    # <-- Usually the most important
    AB: τ = [0.3333, 0.3333]  # Standard high-symmetry stacking
    BA: τ = [0.6667, 0.6667]
    Broken symmetries: Mx, S3, S3^5...

  Pair 2:                    # <-- Alternative configuration
    AB: τ = [0.1667, 0.1667]  # Different stacking pattern
    BA: τ = [0.8333, 0.8333]
    Broken symmetries: Mx, S3, S3^5...
```

## Output Formats

### CIF Format (Default)

By default, SYMPOL2D generates CIF files for each stacking configuration found. These files can be directly visualized in software like VESTA, Jmol, or Materials Studio.

```bash
# Default behavior - creates MoS2_stackings/ directory with CIF files
sympol2d search --formula MoS2 --auto-select

# Output structure:
# MoS2_stackings/
#   ├── MoS2_AA.cif      # Non-polar AA stacking
#   ├── MoS2_AB.cif      # Polar AB stacking
#   └── MoS2_BA.cif      # Polar BA stacking
```

### JSON Format

For programmatic analysis, use JSON output:

```bash
sympol2d search --formula MoS2 --auto-select --output-format json --output results.json
```

## Future Features (Not Yet Implemented)

- **Polarization calculation**: Compute polarization magnitudes
- **Energy estimation**: Estimate relative stabilities of different stackings
- **Batch processing**: Process multiple materials in one run

## Overview

SYMPOL2D is a tool for identifying polar (ferroelectric/antiferroelectric) stacking configurations in 2D bilayer materials through systematic symmetry analysis. It uses a grid-based approach to explore all possible stacking configurations and identifies those that break inversion symmetry.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/sympol2d.git
cd sympol2d

# Install the package
pip install -e .
```

## Database Setup

SYMPOL2D requires the C2DB (Computational 2D Materials Database) to function:

1. **For the included example**: The database is already provided in `example/c2db.db`
2. **For full functionality**: Download the complete C2DB database from [https://cmr.fysik.dtu.dk/c2db/c2db.html](https://cmr.fysik.dtu.dk/c2db/c2db.html)

### Using the Example Database

```bash
# Run with the example database
sympol2d --database example/c2db.db search --uid 1WS2-1

# Or use the runner script
python3 run_sympol2d.py --database example/c2db.db search --uid 1WS2-1
```

## How It Works

### 1. Grid-Based Stacking Exploration

The tool systematically explores stacking configurations using a uniform 2D grid:
- Creates a grid of stacking vectors τ = [τx, τy] in the unit cell
- Each grid point represents a relative shift between layers
- Default: 50×50 grid (2,500 configurations tested)

### 2. Symmetry Analysis

For each stacking configuration:
- Tests which symmetry operations are preserved/broken
- Identifies configurations that break inversion symmetry (polar)
- Determines polarization direction (x, y, z, xy, or general)

### 3. Stacking Classification

- **AA stacking**: τ = [0, 0] - layers perfectly aligned (non-polar)
- **AB/BA stacking**: Related by inversion (BA = 1 - AB)
- Common high-symmetry stackings: τ = [1/3, 1/3] (AB) and [2/3, 2/3] (BA)

## Usage Examples

### Example 1: Finding polar stackings in MoS2

```bash
$ sympol2d search --formula MoS2 --auto-select

Analyzing: MoS2 (MoS2-165798ab5e18)
Layer group: P-3m1
Number of atoms: 3
Auto-estimated interlayer distance: 3.10 Å

Scanning 50x50 grid for stacking configurations...

============================================================
STACKING CONFIGURATIONS FOUND:
============================================================

AA stacking (non-polar):
  τ = [0.0000, 0.0000], d = 3.10 Å
  Preserved symmetries: E, C3, C3^2, Mx, My, Mxy, Mxy-, S3, S3^5, C2

X-POLARIZED PAIRS (12 found):

  Pair 1:
    AB: τ = [0.3333, 0.3333]
    BA: τ = [0.6667, 0.6667]
    Broken symmetries: Mx, S3, S3^5...

  ... and 11 more x-polarized pairs
```

### Example 2: Exploring different materials

```bash
# List available TMDCs
$ sympol2d list --formula WS2
Found 20 materials:
------------------------------------------------------------
UID                            Formula         Layer Group    
------------------------------------------------------------
WS2-4f8d1a2b3c                WS2             P-3m1          
WS2-7a9e5d6f1e                WS2             P-6m2          
...

# Search specific material
$ sympol2d search --uid WS2-4f8d1a2b3c --polar-direction z
```

### Example 3: High-resolution grid search

```bash
# Use finer grid for more accurate results
$ sympol2d search --formula MoS2 --auto-select --grid 100 --output MoS2_fine_grid.json
```

## Output Format

The tool provides:
1. **Console output**: Human-readable summary of findings
2. **JSON output** (optional): Detailed results for further analysis

JSON structure:
```json
{
  "material": {
    "uid": "MoS2-165798ab5e18",
    "formula": "MoS2",
    "layer_group": "P-3m1"
  },
  "grid_size": 50,
  "stackings": {
    "AA": {
      "tau": [0.0, 0.0],
      "preserved_symmetries": ["E", "C3", ...],
      "broken_symmetries": []
    },
    "AB": {
      "tau": [0.3333, 0.3333],
      "preserved_symmetries": ["E", "C3", ...],
      "broken_symmetries": ["Mx", "S3", ...]
    },
    "BA": {
      "tau": [0.6667, 0.6667],
      "preserved_symmetries": ["E", "C3", ...],
      "broken_symmetries": ["Mx", "S3", ...]
    }
  }
}
```

## Polarization Directions

- **x-polar**: Polarization along x-axis (My preserved, Mx broken)
- **y-polar**: Polarization along y-axis (Mx preserved, My broken)
- **z-polar**: Out-of-plane polarization (both Mx and My broken, C2 preserved)
- **xy-polar**: In-plane diagonal polarization
- **general**: Complex polarization pattern

## Technical Details

### Symmetry Operations Tested

The tool tests preservation of:
- Rotations: C2, C3, C4, C6
- Mirrors: Mx, My, Mxy, Mxy-
- Roto-inversions: S3, S6
- Identity: E

### High-Symmetry Stackings

The tool prioritizes common high-symmetry positions:
- 1/6, 1/4, 1/3, 1/2, 2/3, 3/4, 5/6
- Standard TMDC stackings: (1/3, 1/3) and (2/3, 2/3)

### Interlayer Distance Estimation

Auto-estimates based on atomic composition:
- TMDCs (Mo, W): ~3.1 Å
- Group III-VI: ~2.8 Å
- Default: 3.1 Å

## Requirements

- Python 3.7+
- NumPy
- ASE (Atomic Simulation Environment)
- Access to C2DB database file

## License

MIT License

## Citation

If you use SYMPOL2D in your research, please cite:
```
[Your publication details here]
```

## Acknowledgments

This tool uses the Computational 2D Materials Database (C2DB). If you use data from C2DB, please cite:
```
Morten Niklas Gjerding et al 2021 2D Mater. 8 044002
```