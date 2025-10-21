# 2dSYMPOL

**SYMmetry-based prediction of POLarity in 2D bilayers**

A Python tool for identifying sliding ferroelectric AB/BA stacking pairs in 2D bilayer materials using symmetry-based analysis.

## Installation

```bash
# Clone the repository
git clone https://github.com/amal-kishore/SymPol2d.git
cd SymPol2d

# Install dependencies
pip install numpy

# Obtain c2db.db separately from c2db maintainers (requires permission)
# Place the database file in raw/ directory
mkdir -p raw
# Copy your c2db.db file to raw/
```

**Note**: `sqlite3` is usually included with Python.

## Usage

**Input command:**
```bash
python3 run_sympol2d.py search --uid 1MoS2-1 --grid 30 --database raw/c2db.db --export --out-prefix mos2
```

This command executes the SymPol2D search for sliding ferroelectric AB/BA stacking pairs in monolayer MoS₂.

**Arguments:**
- `--uid` specifies the monolayer identifier from the C2DB database
- `--grid 30` performs a 30×30 stacking scan
- `--database` provides the structural data source
- `--export` and `--out-prefix` generate output CIF files for the identified AB and BA stackings

**Output files:** `mos2_AB.cif` and `mos2_BA.cif`

## Methodology

Tests symmetry preservation: `(E + R)τ = n`

For sliding ferroelectrics, 2dSYMPOL finds AB/BA pairs where:
- **AB stacking**: τ ≈ [0.333, 0.333] (e.g., for TMDCs)
- **BA stacking**: τ ≈ [0.667, 0.667] (related by inversion)
- AB and BA have opposite out-of-plane polarization (Pz)

## Documentation

See `MANUAL.md` for detailed usage instructions.

## Acknowledgments

This tool uses the Computational 2D Materials Database (c2db). Please obtain the database separately with permission and cite:
```
Haastrup, S., Strange, M., Pandey, M. et al.
The Computational 2D Materials Database: high-throughput modeling and discovery of atomically thin crystals.
2D Mater. 5, 042002 (2018).
```

## Author

**Amal Kishore**
