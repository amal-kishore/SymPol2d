# 2dSYMPOL Manual

**SYMmetry-based prediction of POLarity in 2D bilayers**

## Table of Contents

1. [Introduction](#introduction)
2. [Installation and Setup](#installation-and-setup)
3. [Command Line Interface](#command-line-interface)
4. [Methodology](#methodology)
5. [Examples](#examples)
6. [Understanding Results](#understanding-results)
7. [Advanced Usage](#advanced-usage)
8. [Troubleshooting](#troubleshooting)
9. [Technical Details](#technical-details)

## Introduction

2dSYMPOL is a computational tool designed to predict polar and non-polar stacking configurations in 2D bilayer systems using crystallographic symmetry analysis. Unlike traditional approaches that rely on expensive density functional theory (DFT) calculations, 2dSYMPOL uses pure symmetry-based matrix operations to rapidly identify promising stacking configurations.

### Key Capabilities

- **Rapid screening**: Analyze thousands of stacking configurations in seconds
- **No DFT required**: Pure symmetry-based approach eliminates computational overhead
- **High-symmetry prioritization**: Automatically identifies physically meaningful stackings
- **Comprehensive analysis**: Find all polar directions (x, y, z, diagonal)
- **Database integration**: Direct access to 16,905 materials from c2db
- **Structure generation**: Create CIF files for visualization and further analysis

### Scientific Background

The polarization in 2D bilayers arises from breaking of inversion symmetry during stacking. Traditional approaches require expensive electronic structure calculations to determine the ground state polarization. 2dSYMPOL bypasses this by using the mathematical relationship between crystallographic symmetry operations and polarization directions.

## Installation and Setup

### Prerequisites

- Python 3.7 or higher
- NumPy
- SQLite3 (usually included with Python)

### Installation Steps

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd slipmat
   ```

2. **Verify Python dependencies**:
   ```bash
   python3 -c "import numpy, sqlite3; print('Dependencies OK')"
   ```

3. **Obtain c2db database**:
   - Download `c2db.db` file separately from c2db maintainers
   - See acknowledgments section for proper citation
   - Place in the project root directory
   - Verify: `ls -la c2db.db` should show the database file

4. **Test installation**:
   ```bash
   python3 run_sympol2d.py --help
   ```

### Directory Structure

```
slipmat/
├── run_sympol2d.py         # Main executable script
├── sympol2d/               # Core package
│   ├── __init__.py
│   ├── cli.py              # Command-line interface
│   ├── symmetry.py         # Symmetry operations
│   ├── scanner.py          # Grid scanning logic
│   ├── c2db_interface.py   # Database interface
│   ├── builder.py          # Bilayer structure builder
│   ├── cif_writer.py       # CIF file generator
│   └── poscar_io.py        # POSCAR file I/O
├── raw/
│   └── c2db.db             # c2db materials database
├── example/                # Example calculations
│   ├── MoS2/              # MoS2 example
│   ├── ReS2/              # ReS2 example
│   ├── hBN/               # hBN example
│   └── WS2/               # WS2 example
├── README.md              # Quick start guide
└── MANUAL.md              # This file
```

## Command Line Interface

2dSYMPOL provides the `search` command to find polar AB/BA stacking pairs for sliding ferroelectrics.

### Search Command

Find the best AB/BA stacking pair for out-of-plane sliding ferroelectricity:

```bash
python3 run_sympol2d.py search [options]
```

**Material selection** (choose one):
- `--uid UID`: c2db material identifier (e.g., `1MoS2-1`)
- `--layer-group GROUP`: Layer group symmetry (e.g., `p-6m2`, `p6mm`, `p2mm`, `pman`)

**Grid scanning**:
- `--grid SIZE`: Grid resolution for tau-space scanning (default: 60)

**Structure export**:
- `--export`: Export bilayer structures (requires `--uid` for CIF or `--poscar` for VASP)
- `--format {cif,poscar}`: Export format (default: `cif`)
- `--out-prefix PREFIX`: Output file prefix (default: `sympol2d`)

**Structure parameters**:
- `--gap DIST`: Interlayer gap in Angstroms (default: 3.1)
- `--poscar FILE`: Monolayer POSCAR file (for VASP format export)
- `--database PATH`: Path to c2db database (default: `raw/c2db.db`)

**Advanced options**:
- `--allow-nonflipping`: Export even if AB/BA don't have opposite Pz

### Examples of Common Usage

```bash
# Basic search with CIF export from c2db
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 30 --database raw/c2db.db --export

# Custom interlayer gap
python3 run_sympol2d.py search --uid 1MoS2-1 --gap 3.2 --export --database raw/c2db.db

# Export with custom output prefix
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 60 --export --out-prefix mos2/mos2

# Using monolayer POSCAR for VASP format
python3 run_sympol2d.py search --layer-group p-6m2 --poscar POSCAR_mono --gap 3.1 --export --format poscar

# Rectangular system (no Pz flip expected)
python3 run_sympol2d.py search --layer-group pman --poscar BP_mono.vasp --export --allow-nonflipping --format poscar

# Quick search without export
python3 run_sympol2d.py search --uid 1WS2-1 --grid 30 --database raw/c2db.db
```

## Methodology

### Mathematical Framework

2dSYMPOL is based on the symmetry preservation condition for stacking vectors:

```
(E + R)τ = n
```

Where:
- **E**: 2×2 identity matrix
- **R**: 2×2 matrix representation of a symmetry operation
- **τ**: 2D stacking vector [τₓ, τᵧ]
- **n**: Integer lattice vector [nₓ, nᵧ]

**Interpretation**:
- If this equation has integer solutions, symmetry operation R is preserved
- If no integer solution exists, the symmetry is broken
- Different combinations of preserved/broken symmetries determine polarization direction

### Classification Logic

1. **AA stacking (non-polar)**:
   - All symmetries preserved
   - τ typically close to [0, 0] or high-symmetry positions

2. **AB stacking (polar)**:
   - Specific symmetries broken (mirrors, inversions)
   - Classified by which symmetries remain intact

3. **BA stacking (polar)**:
   - Partner to AB with opposite polarization
   - Related by inversion: BA = 1 - AB (modulo lattice)

### Polarization Direction Assignment

- **x-polar**: My mirror preserved, Mx mirror broken
- **y-polar**: Mx mirror preserved, My mirror broken  
- **z-polar**: Both Mx and My broken, but C₂ rotation preserved
- **xy-polar**: Diagonal mirrors (Mxy) broken
- **general**: Complex symmetry breaking pattern

### Grid Scanning Strategy

1. **Dense grid generation**: Creates N×N grid over unit cell (0,0) to (1,1)
2. **Symmetry testing**: Each grid point tested against all layer group operations
3. **High-symmetry prioritization**: Points at simple fractions (1/6, 1/4, 1/3, 1/2, 2/3, 3/4, 5/6) ranked higher
4. **AB-BA pairing**: Inverted stackings identified and paired
5. **Physical validation**: Ensures AB/BA pairs represent true opposite polarities

## Examples

### Example 1: MoS2 Analysis with CIF Export

MoS2 (molybdenum disulfide) is a prototypical transition metal dichalcogenide with layer group p-6m2.

```bash
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 30 --database raw/c2db.db --export
```

**Expected output**:
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

Extracting structure from c2db database...
Loaded: MoS2 (p-6m2), 3 atoms

Exported bilayer CIF structures:
  AB: sympol2d_AB.cif
  BA: sympol2d_BA.cif
```

**Analysis**:
- MoS2 shows out-of-plane polarization (sliding ferroelectric)
- Standard AB/BA stackings at τ = [1/3, 1/3] and [2/3, 2/3]
- These correspond to standard 2H-type stacking in TMDCs
- CIF files can be visualized in VESTA, Materials Studio, or similar tools

### Example 2: Custom Output Location

```bash
python3 run_sympol2d.py search --uid 1WS2-1 --grid 60 --export --out-prefix example/WS2/ws2 --database raw/c2db.db
```

This will create `example/WS2/ws2_AB.cif` and `example/WS2/ws2_BA.cif`.

### Example 3: POSCAR Format Export

For VASP calculations, you can export in POSCAR format if you have the monolayer structure:

```bash
python3 run_sympol2d.py search --layer-group p-6m2 --poscar POSCAR_mono --gap 3.1 --export --format poscar --out-prefix bilayer
```

This will create `bilayer_AB.vasp` and `bilayer_BA.vasp`.

## Understanding Results

### Output Interpretation

1. **Grid Scanning Results**:
   - Number of grid points: N×N where N is your `--grid` value
   - Number of stackings found: Configurations that break inversion symmetry

2. **Best AB/BA Pair**:
   - **Layer group**: Crystallographic symmetry of the monolayer
   - **Formula**: Chemical composition
   - **AB stacking**: τ vector for AB configuration
   - **BA stacking**: τ vector for BA configuration (typically related to AB by inversion)

3. **Polarization Check**:
   - **✓ AB and BA have opposite Pz**: Material is a sliding ferroelectric
   - **✗ AB and BA do NOT have opposite Pz**: Not a sliding ferroelectric (some layer groups)

4. **Exported Files**:
   - CIF format: `{prefix}_AB.cif` and `{prefix}_BA.cif`
   - POSCAR format: `{prefix}_AB.vasp` and `{prefix}_BA.vasp`

### Quality Indicators

**High-quality results**:
- AB/BA pairs at simple fractions (1/3, 1/2, 2/3)
- Clear message about Pz flipping
- Reasonable number of candidate stackings (10-1000)

**Potential issues**:
- τ values very close to grid artifacts (e.g., 0.02, 0.98) - try higher grid resolution
- Very few stackings (<5) - may need higher grid resolution
- Thousands of stackings (>2000) - may indicate numerical issues

### Physical Interpretation

**AB vs BA stackings**:
- AB: One specific lateral stacking arrangement
- BA: Related by inversion symmetry: τ_BA ≈ 1 - τ_AB (modulo 1)
- For TMDCs: Often corresponds to "chalcogen over metal" vs "metal over chalcogen"
- Sliding between AB and BA flips the out-of-plane polarization

**Out-of-plane polarization (Pz)**:
- Most common in van der Waals materials with broken inversion symmetry
- AB and BA states have opposite Pz → sliding ferroelectric
- Useful for non-volatile memory, sensors, and energy harvesting devices

## Advanced Usage

### Custom Grid Scanning

For high-precision work, increase grid density:

```bash
# Ultra-high resolution (10,000 points)
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 100 --database raw/c2db.db

# Quick screening (100 points)
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 10 --database raw/c2db.db

# Standard (default: 60×60 = 3600 points)
python3 run_sympol2d.py search --uid 1MoS2-1 --database raw/c2db.db
```

**Grid density guidelines**:
- Grid 10-20: Quick screening, may miss fine features
- Grid 30-60: Standard analysis, good balance (default: 60)
- Grid 70-100: High precision, computational cost increases

### Interlayer Gap Effects

```bash
# Standard van der Waals gap
python3 run_sympol2d.py search --uid 1WS2-1 --gap 3.1 --export --database raw/c2db.db

# Compressed bilayers
python3 run_sympol2d.py search --uid 1WS2-1 --gap 2.5 --export --database raw/c2db.db

# Expanded interlayers
python3 run_sympol2d.py search --uid 1WS2-1 --gap 4.0 --export --database raw/c2db.db
```

Note: The `--gap` parameter only affects the exported structure geometry, not the symmetry-based analysis.

### Batch Processing

For systematic studies, create bash scripts:

```bash
#!/bin/bash
# batch_analysis.sh

materials=("1MoS2-1" "1WS2-1" "1WSe2-1")
for mat in "${materials[@]}"; do
    echo "Processing $mat..."
    python3 run_sympol2d.py search --uid "$mat" --grid 30 \
        --export --out-prefix "results/${mat}/${mat}" \
        --database raw/c2db.db
done
```

This will create directory structure:
```
results/
├── 1MoS2-1/
│   ├── 1MoS2-1_AB.cif
│   └── 1MoS2-1_BA.cif
├── 1WS2-1/
│   ├── 1WS2-1_AB.cif
│   └── 1WS2-1_BA.cif
└── 1WSe2-1/
    ├── 1WSe2-1_AB.cif
    └── 1WSe2-1_BA.cif
```

## Troubleshooting

### Common Issues

**1. "Database file not found"**
```bash
# Verify database location (default path)
ls -la raw/c2db.db

# Or specify custom path
python3 run_sympol2d.py search --uid 1MoS2-1 --database /path/to/c2db.db

# Check file permissions
chmod 644 raw/c2db.db
```

**2. "No suitable AB/BA pair found"**
- Some layer groups may not have inversion-breaking stackings
- Try increasing grid resolution: `--grid 80`
- Check if the layer group is compatible with sliding ferroelectricity

**3. "CIF export requires --uid"**
- CIF export needs the c2db database to extract atomic structure
- Either provide `--uid` with `--database`, or use POSCAR format instead:
```bash
python3 run_sympol2d.py search --layer-group p-6m2 --poscar POSCAR_mono --format poscar --export
```

**4. "POSCAR file not found"**
```bash
# Verify POSCAR exists
ls -la POSCAR_mono

# Use absolute or relative path
python3 run_sympol2d.py search --layer-group p-6m2 --poscar /full/path/to/POSCAR_mono --export --format poscar
```

**5. ImportError or module issues**
```bash
# Verify Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Check dependencies
python3 -c "import numpy; print('NumPy OK')"

# Or run directly as module
python3 -m sympol2d.cli search --uid 1MoS2-1 --database raw/c2db.db
```

**6. "provide --layer-group or a valid --uid"**
- You must specify either `--uid` (with database) OR `--layer-group`
- With `--uid`, layer group is automatically determined from database
- Without database, you must manually specify `--layer-group`

### Performance Optimization

**Memory usage**:
- Large grids (>80) require more RAM
- Consider smaller grids (30-60) for batch processing

**Speed optimization**:
- Smaller grids complete faster (grid 30 vs grid 100)
- CIF export is fast (direct from database)
- POSCAR export requires monolayer file reading

### Getting Help

**Debug mode**:
Add print statements in source code for detailed analysis:

```python
# In scanner.py, add debugging
print(f"Testing tau = {tau}, preserved = {preserved_ops}")
```

**Verify symmetry operations**:
```python
# Check layer group operations
from sympol2d.symmetry import LayerGroupSymmetry
symmetry = LayerGroupSymmetry('p-6m2')
print(symmetry.operations)
```

## Technical Details

### Supported Layer Groups

Currently implemented:
- **p-6m2**: Hexagonal (TMDCs like MoS2, WS2)
- **p-4m2**: Square lattice
- **p-3m1**: Triangular lattice
- **p-mmm**: Rectangular lattice

Adding new layer groups requires:
1. Define symmetry operations in `sympol2d/symmetry.py`
2. Add to `LAYER_GROUP_OPERATIONS` dictionary
3. Test with known materials

### Matrix Representations

Symmetry operations use 2×2 matrices:
```python
OPERATIONS = {
    'E':    np.array([[1,  0], [0,  1]]),   # Identity
    'C2':   np.array([[-1, 0], [0, -1]]),   # 180° rotation
    'Mx':   np.array([[1,  0], [0, -1]]),   # Mirror x
    'My':   np.array([[-1, 0], [0,  1]]),   # Mirror y
    'C6':   np.array([[1/2, -sqrt(3)/2], [sqrt(3)/2, 1/2]]),  # 60° rotation
    ...
}
```

### Algorithm Complexity

- **Grid scanning**: O(N²) where N is grid size
- **Symmetry testing**: O(S) where S is number of symmetry operations
- **Pair finding**: O(P²) where P is number of polar stackings
- **Overall**: O(N² × S + P²)

For typical materials:
- N = 50, S = 8, P = 100-1000
- Runtime: 1-10 seconds on modern hardware

### Numerical Precision

**Grid spacing effects**:
- Grid 50: spacing = 0.02 (2% of unit cell)
- High-symmetry fractions may be missed if spacing too large
- Tolerance for matching: typically 2× grid spacing

**Floating point considerations**:
- Use `np.allclose()` with appropriate tolerance
- Modular arithmetic for stacking vectors
- Integer solutions checked with small tolerance (1e-10)

### Database Schema

c2db database contains:
- **uid**: Unique material identifier
- **layer_group**: Crystallographic symmetry
- **natoms**: Number of atoms in unit cell
- **numbers**: Atomic numbers array
- Additional structural and electronic properties

### Extension Points

**New polarization types**:
1. Define new classification logic in `_determine_polar_direction()`
2. Add to CLI choices in `--polar-direction`
3. Update documentation

**Alternative databases**:
1. Implement new interface inheriting from base class
2. Ensure consistent Material2D dataclass format
3. Update CLI to accept database type parameter

**Export formats**:
1. Add new output formats in `search_material()`
2. Implement structure writers (VASP, Quantum ESPRESSO, etc.)
3. Consider integration with ASE (Atomic Simulation Environment)

## Acknowledgments

This tool uses the Computational 2D Materials Database (c2db). Please obtain the database separately and cite:

```
Haastrup, S., Strange, M., Pandey, M. et al. 
The Computational 2D Materials Database: high-throughput modeling and discovery of atomically thin crystals. 
2D Mater. 5, 042002 (2018).
```

## Author

**Amal Kishore**

---

*For additional questions or issues, please consult the source code or contact the author.*