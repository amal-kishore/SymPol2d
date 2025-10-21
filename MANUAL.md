# 2dSYMPOL Manual

**SYMmetry-based prediction of POLarity in 2D bilayers**

## Introduction

2dSYMPOL finds polar AB/BA stacking pairs in 2D bilayer materials using symmetry analysis. It identifies sliding ferroelectric configurations where AB and BA stackings have opposite out-of-plane polarization.

## Installation

**Prerequisites**: Python 3.7+, NumPy, SQLite3

**Setup**:
```bash
git clone <repository-url>
cd slipmat
python3 -c "import numpy, sqlite3; print('OK')"
```

**Database** (required for CIF export):
- Download `c2db.db` separately from c2db maintainers (requires permission)
- Place in `raw/` directory

## Usage

### Basic Command

```bash
python3 run_sympol2d.py search --uid <UID> --database raw/c2db.db [--export]
```

### Key Arguments

**Material selection** (choose one):
- `--uid UID`: c2db material ID (e.g., `1MoS2-1`)
- `--layer-group GROUP`: Layer group if no database (e.g., `p-6m2`)

**Export options**:
- `--export`: Export bilayer structures
- `--format {cif,poscar}`: Output format (default: `cif`)
- `--out-prefix PREFIX`: File prefix (default: `sympol2d`)
- `--gap DIST`: Interlayer gap in Å (default: 3.1)

**Scanning**:
- `--grid SIZE`: Resolution (default: 60, range: 10-100)
- `--database PATH`: Database path (default: `raw/c2db.db`)

### Examples

```bash
# Export MoS2 bilayer CIF files
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 30 --database raw/c2db.db --export

# Custom output location
python3 run_sympol2d.py search --uid 1WS2-1 --export --out-prefix results/ws2 --database raw/c2db.db

# POSCAR format (requires monolayer file)
python3 run_sympol2d.py search --layer-group p-6m2 --poscar POSCAR_mono --export --format poscar

# High-resolution scan
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 100 --database raw/c2db.db
```

## Output

**Console output**:
```
Scanning 30×30 grid (900 points) for layer group 'p-6m2'...
  → Found 126 stackings with broken inversion symmetry

============================================================
BEST AB/BA STACKING PAIR
============================================================
Layer group: p-6m2
Formula: MoS2

AB stacking: τ = [0.333333, 0.333333]
BA stacking: τ = [0.666667, 0.666667]

Out-of-plane polarization (Pz):
  ✓ AB and BA have opposite Pz (sliding ferroelectric)
============================================================

Exported bilayer CIF structures:
  AB: sympol2d_AB.cif
  BA: sympol2d_BA.cif
```

**Exported files**:
- CIF: `{prefix}_AB.cif`, `{prefix}_BA.cif`
- VASP: `{prefix}_AB.vasp`, `{prefix}_BA.vasp`

## Methodology

Tests symmetry preservation: `(E + R)τ = n`
- **E**: Identity matrix
- **R**: Symmetry operation matrix
- **τ**: Stacking vector [τₓ, τᵧ]
- **n**: Integer lattice vector

If no integer solution exists, symmetry is broken. The tool:
1. Scans a grid of stacking vectors
2. Tests which symmetries are preserved/broken
3. Identifies AB/BA pairs related by inversion
4. Prioritizes high-symmetry configurations (τ = 1/3, 1/2, 2/3, etc.)

For sliding ferroelectrics, AB and BA must have opposite out-of-plane polarization.

## Troubleshooting

**Database not found**:
```bash
ls -la raw/c2db.db
# Place database in raw/ directory
```

**No AB/BA pair found**:
- Try higher grid resolution: `--grid 80`
- Some materials don't support sliding ferroelectricity

**CIF export requires --uid**:
- CIF needs database for atomic structure
- Alternative: use POSCAR format with `--poscar` and `--format poscar`

**Module import errors**:
```bash
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
python3 -c "import numpy; print('NumPy OK')"
```

## Grid Resolution Guidelines

- **10-20**: Quick screening (may miss features)
- **30-60**: Standard (default: 60)
- **70-100**: High precision (slower)

## Supported Layer Groups

- **p-6m2**: Hexagonal (MoS2, WS2)
- **p-4m2**: Square lattice
- **p-3m1**: Triangular lattice
- **p-mmm, pman**: Rectangular lattice
- **p-1**: Triclinic

## Acknowledgments

This tool uses the Computational 2D Materials Database (c2db).
**IMPORTANT**: Obtain c2db.db separately with permission. Cite:

```
Haastrup, S., Strange, M., Pandey, M. et al.
The Computational 2D Materials Database: high-throughput modeling
and discovery of atomically thin crystals.
2D Mater. 5, 042002 (2018).
```

**Author**: Amal Kishore
