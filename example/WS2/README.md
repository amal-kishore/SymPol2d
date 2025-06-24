# WS2 Bilayer Analysis Example

This directory contains a comprehensive example of using 2dSYMPOL to analyze WS2 bilayer stacking configurations.

## Material Information

- **Material**: WS2 (Tungsten Disulfide)
- **c2db UID**: 1WS2-1
- **Layer Group**: p-6m2 (hexagonal symmetry)
- **Monolayer Structure**: 3 atoms (1 W + 2 S)
- **Lattice Parameters**: a = b = 3.186 Å, c = 18.146 Å, γ = 120°

## Analysis Results

### Symmetry Operations
WS2 has p-6m2 layer group symmetry with operations:
- E (identity)
- C6, C3, C2 (rotations: 60°, 120°, 180°)
- C3², C6⁵ (additional rotations)
- Mx, My, Mxy (mirror planes)

### Stacking Configurations Found

#### Grid Scan (30×30)
- **AA stackings**: 1 configuration (non-polar)
- **X-polarized pairs**: 28 AB-BA pairs
- **Y-polarized pairs**: 28 AB-BA pairs  
- **Z-polarized pairs**: 392 AB-BA pairs

#### Representative Stackings
1. **AA (Non-polar)**: τ = [0.0000, 0.0000]
   - All symmetries preserved
   - Interlayer distance: 3.10 Å

2. **AB (Z-polar)**: τ = [0.0333, 0.0333]
   - Broken symmetries: E, C6, C3, C2, C3², C6⁵, Mx, My, Mxy
   - Out-of-plane polarization

3. **BA (Z-polar opposite)**: τ = [0.9667, 0.9667]
   - Same broken symmetries as AB
   - Opposite out-of-plane polarization

## Files Generated

### Analysis Results
- `WS2_analysis_full.txt`: Complete analysis output
- `WS2_all_stackings.json`: JSON data of all configurations

### CIF Files (in withcif/ directory)
- `1WS2-1.cif`: Original monolayer structure
- `WS2_bilayer_AA.cif`: AA stacking bilayer
- `WS2_bilayer_AB.cif`: AB stacking bilayer (z-polar)
- `WS2_bilayer_BA.cif`: BA stacking bilayer (z-polar opposite)
- `generate_bilayer_cifs.py`: Script to generate bilayer structures

### Bilayer Structure Details
- **Unit Cell**: c = 21.246 Å (monolayer + 3.1 Å gap)
- **Total Atoms**: 6 (2 layers × 3 atoms/layer)
- **Interlayer Distance**: 3.1 Å (typical vdW distance)

## Key Findings

1. **Polarity**: WS2 bilayers show strong out-of-plane (z) polarization capability
2. **Abundance**: 392 different z-polarized AB-BA pairs possible
3. **Symmetry Breaking**: Polar stackings break multiple symmetries while preserving some rotational order
4. **Structure**: AB and BA pairs have τ vectors related by inversion symmetry

## Usage Examples

### Basic Analysis
```bash
python3 run_2dsympol.py search --uid 1WS2-1 --grid 30
```

### Filter by Polarization
```bash
python3 run_2dsympol.py search --uid 1WS2-1 --polar-direction z
```

### Custom Interlayer Distance
```bash
python3 run_2dsympol.py search --uid 1WS2-1 --interlayer-distance 3.5
```

## Physical Interpretation

The z-polarized AB and BA stackings represent bilayer configurations where:
- **AB**: Top layer shifted by (+0.0333, +0.0333) relative to bottom
- **BA**: Top layer shifted by (-0.0333, -0.0333) relative to bottom (equivalent to +0.9667, +0.9667 due to periodicity)

These configurations break inversion symmetry, leading to:
- Electric dipole moment perpendicular to layers
- Potential for ferroelectric switching between AB and BA
- Different electronic properties compared to AA stacking

## Applications

These bilayer structures can be used for:
1. **DFT Calculations**: Input structures for electronic property calculations
2. **Experimental Validation**: Comparison with STM/AFM stacking measurements  
3. **Device Design**: Understanding polar bilayer behavior in van der Waals heterostructures
4. **Machine Learning**: Training data for stacking property prediction models