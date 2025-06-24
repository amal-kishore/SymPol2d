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
   python3 run_2dsympol.py --help
   ```

### Directory Structure

```
slipmat/
├── run_2dsympol.py          # Main executable script
├── sympol2d/               # Core package
│   ├── __init__.py
│   ├── cli.py              # Command-line interface
│   ├── symmetry.py         # Symmetry operations
│   ├── scanner.py          # Grid scanning logic
│   ├── c2db_interface.py   # Database interface
│   └── utils.py            # Utility functions
├── c2db.db                 # c2db materials database
├── example/                # Example calculations
│   └── WS2/               # WS2 example with CIF files
├── README.md              # Quick start guide
└── MANUAL.md              # This file
```

## Command Line Interface

2dSYMPOL provides two main commands: `search` and `list`.

### Search Command

Find polar stackings for a specific material:

```bash
python3 run_2dsympol.py search [options]
```

**Required arguments** (one of):
- `--uid UID`: c2db material identifier (e.g., `1WS2-1`)
- `--formula FORMULA`: Chemical formula (e.g., `MoS2`)

**Optional arguments**:
- `--grid SIZE`: Grid density (default: 50, range: 10-100)
- `--polar-direction DIR`: Filter by direction (`x`, `y`, `z`, `xy`, `general`, `all`)
- `--interlayer-distance DIST`: Interlayer separation in Å (default: auto-estimate)
- `--output FILE`: Save results to JSON file
- `--auto-select`: Auto-select first match when multiple materials found

### List Command

Browse available materials:

```bash
python3 run_2dsympol.py list [options]
```

**Options**:
- `--formula FORMULA`: Filter by chemical formula
- `--layer-group GROUP`: Filter by layer group symmetry
- `--layer-groups`: Show all available layer groups
- `--limit N`: Maximum results to display (default: 20)

### Examples of Common Usage

```bash
# Basic search by material ID
python3 run_2dsympol.py search --uid 1WS2-1

# Search by formula with auto-selection
python3 run_2dsympol.py search --formula MoS2 --auto-select

# High-resolution scan for z-polar stackings
python3 run_2dsympol.py search --uid 1WS2-1 --grid 80 --polar-direction z

# Custom interlayer distance
python3 run_2dsympol.py search --uid 1WS2-1 --interlayer-distance 3.5

# Save detailed results
python3 run_2dsympol.py search --uid 1WS2-1 --output ws2_results.json

# Browse materials by composition
python3 run_2dsympol.py list --formula WS2

# Show all available layer groups
python3 run_2dsympol.py list --layer-groups
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

### Example 1: WS2 Analysis

WS2 (tungsten disulfide) is a prototypical transition metal dichalcogenide with layer group p-6m2.

```bash
python3 run_2dsympol.py search --uid 1WS2-1 --grid 50
```

**Expected output**:
```
Analyzing: WS2 (1WS2-1)
Layer group: p-6m2
Number of atoms: 3
Auto-estimated interlayer distance: 3.10 Å

Scanning 50x50 grid for stacking configurations...

============================================================
STACKING CONFIGURATIONS FOUND:
============================================================

AA stacking (non-polar):
  τ = [0.0000, 0.0000], d = 3.10 Å
  Preserved symmetries: E, C6, C3, C2, C3^2, C6^5, Mx, My, Mxy

Z-POLARIZED PAIRS (392 found):

  Pair 1:
    AB: τ = [0.3333, 0.3333]
    BA: τ = [0.6667, 0.6667]
    Broken symmetries: Mx, My, Mxy
```

**Analysis**:
- WS2 shows strong z-polarization (out-of-plane)
- Standard AB/BA stackings at τ = [1/3, 1/3] and [2/3, 2/3]
- These correspond to standard 2H-type stacking in TMDCs
- 392 total z-polar pairs found, but high-symmetry ones prioritized

### Example 2: MoS2 Comparison

```bash
python3 run_2dsympol.py search --formula MoS2 --auto-select --polar-direction z
```

MoS2 should show similar behavior to WS2 due to identical layer group symmetry.

### Example 3: Square Lattice Material

```bash
python3 run_2dsympol.py search --formula <square-lattice-material> --polar-direction all
```

Square lattice materials (layer group p-4m2) typically show:
- x-polar and y-polar pairs
- z-polar configurations
- Different symmetry breaking patterns than hexagonal systems

## Understanding Results

### Output Interpretation

1. **Material Information**:
   - UID: Unique identifier in c2db
   - Formula: Chemical composition
   - Layer group: Crystallographic symmetry
   - Atom count: Number of atoms in unit cell

2. **Scanning Parameters**:
   - Grid size: Density of stacking vector sampling
   - Interlayer distance: Vertical separation between layers

3. **AA Stacking**:
   - Always shown first
   - τ ≈ [0, 0] for most materials
   - Lists all preserved symmetries

4. **Polar Pairs**:
   - Grouped by polarization direction
   - Ordered by symmetry importance
   - Shows number of pairs found
   - AB/BA relationship explicitly shown

### Quality Indicators

**High-quality results**:
- AB/BA pairs at simple fractions (1/3, 1/2, 2/3)
- Clear symmetry breaking patterns
- Reasonable number of pairs (not too many/few)

**Potential issues**:
- τ values very close to grid artifacts (e.g., 0.02, 0.98)
- Thousands of polar pairs (may indicate numerical issues)
- No polar pairs found (check layer group compatibility)

### Physical Interpretation

**AB vs BA stackings**:
- AB: Often corresponds to "chalcogen over metal" in TMDCs
- BA: Opposite arrangement with inverted polarization
- Magnitude depends on atomic charges and positions

**Polarization directions**:
- **z-polar**: Most common in van der Waals materials
- **x/y-polar**: Often seen in materials with rectangular unit cells
- **xy-polar**: Diagonal polarization, less common

## Advanced Usage

### Custom Grid Scanning

For high-precision work, increase grid density:

```bash
# Ultra-high resolution (10,000 points)
python3 run_2dsympol.py search --uid 1WS2-1 --grid 100

# Quick screening (100 points)
python3 run_2dsympol.py search --uid 1WS2-1 --grid 10
```

**Grid density guidelines**:
- Grid 10-20: Quick screening, may miss fine features
- Grid 30-50: Standard analysis, good balance
- Grid 60-100: High precision, computational cost increases

### Interlayer Distance Effects

```bash
# van der Waals materials
python3 run_2dsympol.py search --uid 1WS2-1 --interlayer-distance 3.1

# Compressed bilayers
python3 run_2dsympol.py search --uid 1WS2-1 --interlayer-distance 2.5

# Expanded interlayers
python3 run_2dsympol.py search --uid 1WS2-1 --interlayer-distance 4.0
```

Note: Interlayer distance affects the classification for structure generation but not the symmetry-based polar analysis.

### Batch Processing

For systematic studies, create bash scripts:

```bash
#!/bin/bash
# batch_analysis.sh

materials=("1WS2-1" "1MoS2-3" "1WSe2-1")
for mat in "${materials[@]}"; do
    echo "Processing $mat..."
    python3 run_2dsympol.py search --uid "$mat" --output "${mat}_results.json"
done
```

### JSON Output Analysis

Results can be saved in JSON format for further analysis:

```bash
python3 run_2dsympol.py search --uid 1WS2-1 --output results.json
```

JSON structure:
```json
{
  "material": {
    "uid": "1WS2-1",
    "formula": "WS2",
    "layer_group": "p-6m2"
  },
  "grid_size": 50,
  "stackings": {
    "AA": {
      "tau": [0.0, 0.0],
      "preserved_symmetries": ["E", "C6", ...],
      "broken_symmetries": []
    },
    "AB": {
      "tau": [0.3333, 0.3333],
      "preserved_symmetries": ["E", "C6", ...],
      "broken_symmetries": ["Mx", "My", ...]
    }
  }
}
```

## Troubleshooting

### Common Issues

**1. "Material not found in database"**
```bash
# Check if UID exists
python3 run_2dsympol.py list --formula WS2

# Use exact UID from list output
python3 run_2dsympol.py search --uid 1WS2-1
```

**2. "No polar pairs found"**
- Check layer group compatibility
- Some materials may only have AA stackings
- Try increasing grid size: `--grid 80`

**3. "Too many polar pairs (>1000)"**
- May indicate numerical precision issues
- Try smaller grid size: `--grid 30`
- Check if material has unusual symmetry

**4. "Database file not found"**
```bash
# Verify database location
ls -la c2db.db

# Check file permissions
chmod 644 c2db.db
```

**5. ImportError or module issues**
```bash
# Verify Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Check dependencies
python3 -c "import numpy; print('NumPy OK')"
```

### Performance Optimization

**Memory usage**:
- Large grids (>80) require significant RAM
- Consider smaller grids for batch processing

**Speed optimization**:
- Use `--polar-direction z` if only out-of-plane polarization needed
- Smaller grids complete faster
- JSON output adds minimal overhead

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