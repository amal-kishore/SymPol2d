#!/usr/bin/env python3
"""
Analyze black phosphorene (4P-1) for out-of-plane polarization
"""

import sys
import os
sys.path.insert(0, os.path.abspath('../..'))

import numpy as np
from sympol2d.c2db_interface import C2DBInterface
from sympol2d.symmetry import LayerGroupSymmetry
from sympol2d.scanner import StackingScanner
from sympol2d.cif_writer import write_bilayer_cif
from sympol2d.utils import estimate_interlayer_distance

# Connect to database
db = C2DBInterface('../../raw/c2db.db')

# Get black phosphorene (4P-1)
material = db.get_material_by_uid('4P-1')

if not material:
    print("Error: Black phosphorene (4P-1) not found in database")
    sys.exit(1)

print('='*70)
print('BLACK PHOSPHORENE (4P-1) STRUCTURE ANALYSIS')
print('='*70)

print(f'\nFormula: {material.formula}')
print(f'Layer group: {material.layer_group}')
print(f'Number of atoms: {material.natoms}')

print(f'\nLattice parameters:')
a = np.linalg.norm(material.lattice[0])
b = np.linalg.norm(material.lattice[1])
c = np.linalg.norm(material.lattice[2])
print(f'  a = {a:.4f} Å')
print(f'  b = {b:.4f} Å')
print(f'  c = {c:.4f} Å')

print(f'\nLattice vectors:')
for i, vec in enumerate(material.lattice):
    print(f'  a{i+1}: [{vec[0]:.4f}, {vec[1]:.4f}, {vec[2]:.4f}]')

print(f'\nAtomic positions (Cartesian coordinates):')
symbols = material.get_chemical_symbols()
for i, (sym, pos) in enumerate(zip(symbols, material.positions)):
    print(f'  {i+1}. {sym}: ({pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f})')

# Get fractional coordinates
inv_lattice = np.linalg.inv(material.lattice)
frac_positions = material.positions @ inv_lattice.T
print(f'\nAtomic positions (Fractional coordinates):')
for i, (sym, pos) in enumerate(zip(symbols, frac_positions)):
    print(f'  {i+1}. {sym}: ({pos[0]:.6f}, {pos[1]:.6f}, {pos[2]:.6f})')

# Check z-coordinates to understand the layer structure
z_coords = material.positions[:, 2]
unique_z = np.unique(np.round(z_coords, 4))
print(f'\nUnique z-coordinates: {unique_z}')
print(f'Z-thickness of monolayer: {np.max(z_coords) - np.min(z_coords):.4f} Å')

# Analyze layer group symmetry
print(f'\n{"="*70}')
print('SYMMETRY ANALYSIS')
print('='*70)

lg = LayerGroupSymmetry(material.layer_group)
print(f'\nLayer group: {material.layer_group}')
print(f'Symmetry operations present:')
for op in lg.operations:
    print(f'  - {op.name} ({op.type})')

# Check key symmetries
has_inversion = any(op.name == 'C2' for op in lg.operations)
mirror_ops = [op.name for op in lg.operations if 'M' in op.name]

print(f'\nKey symmetry properties:')
print(f'  Has inversion symmetry (C2): {has_inversion}')
print(f'  Has mirror planes: {len(mirror_ops) > 0}')
if mirror_ops:
    print(f'    Mirror operations: {", ".join(mirror_ops)}')

# Now analyze bilayer stackings for polarization
print(f'\n{"="*70}')
print('BILAYER STACKING ANALYSIS FOR OUT-OF-PLANE POLARIZATION')
print('='*70)

# Estimate interlayer distance
d_estimate = estimate_interlayer_distance(material.formula)
print(f'\nEstimated interlayer distance: {d_estimate:.2f} Å')

# Create scanner
scanner = StackingScanner(
    material,
    interlayer_distance=d_estimate,
    grid_size=10  # 10x10 grid for shifts
)

# Scan for polar stackings
print(f'\nScanning {scanner.grid_size}x{scanner.grid_size} = {scanner.grid_size**2} stacking configurations...')
results = scanner.scan_stackings()

# Analyze results
polar_stackings = [r for r in results if r.breaks_inversion]
print(f'\nResults:')
print(f'  Total configurations scanned: {len(results)}')
print(f'  Polar configurations (broken inversion): {len(polar_stackings)}')
print(f'  Non-polar configurations: {len(results) - len(polar_stackings)}')

if polar_stackings:
    print(f'\n{"="*70}')
    print('POLAR STACKINGS FOUND - OUT-OF-PLANE POLARIZATION POSSIBLE!')
    print('='*70)

    print(f'\nTop polar stackings by polarization magnitude:')
    # Sort by polarization magnitude
    polar_stackings.sort(key=lambda x: np.linalg.norm(x.polarization) if x.polarization is not None else 0, reverse=True)

    for i, stacking in enumerate(polar_stackings[:5], 1):
        shift = stacking.shift
        pol = stacking.polarization
        if pol is not None:
            pol_z = pol[2] if len(pol) > 2 else 0
            pol_mag = np.linalg.norm(pol)
            print(f'\n  {i}. Shift: ({shift[0]:.3f}, {shift[1]:.3f})')
            print(f'     Fractional shift: ({shift[0]/a:.3f}, {shift[1]/b:.3f})')
            print(f'     Polarization: [{pol[0]:.4f}, {pol[1]:.4f}, {pol[2]:.4f}] e·Å')
            print(f'     Out-of-plane (z) component: {pol_z:.4f} e·Å')
            print(f'     Total magnitude: {pol_mag:.4f} e·Å')
            print(f'     Broken symmetries: {", ".join(stacking.broken_symmetries)}')

    # Save the best polar stacking
    best_polar = polar_stackings[0]
    shift = best_polar.shift

    # Create shifted top layer
    top_positions = material.positions.copy()
    top_positions[:, 0] += shift[0]
    top_positions[:, 1] += shift[1]
    top_positions[:, 2] += d_estimate

    write_bilayer_cif(
        'best_polar_stacking.cif',
        material.lattice,
        material.lattice,
        material.positions,
        top_positions,
        material.numbers,
        material.numbers,
        interlayer_distance=d_estimate
    )
    print(f'\nBest polar stacking saved to best_polar_stacking.cif')

else:
    print(f'\n{"="*70}')
    print('NO POLAR STACKINGS FOUND')
    print('='*70)

    print('\nWhy black phosphorene does not show out-of-plane polarization:')
    print('\n1. SYMMETRY CONSTRAINTS:')
    print(f'   - Layer group {material.layer_group} has the following symmetry operations:')
    for op in lg.operations:
        print(f'     • {op.name}: {op.type}')

    if has_inversion:
        print('\n2. INVERSION SYMMETRY:')
        print('   - The monolayer has inversion symmetry (C2 rotation)')
        print('   - This constrains the possible stackings')

    if mirror_ops:
        print('\n3. MIRROR PLANES:')
        print(f'   - The structure has mirror symmetry: {", ".join(mirror_ops)}')
        print('   - Mirror planes parallel to the layer suppress out-of-plane polarization')

    print('\n4. STRUCTURAL FACTORS:')
    print('   - Black phosphorene has a puckered honeycomb structure')
    print('   - All P atoms are equivalent by symmetry')
    print('   - The structure maintains its non-polar character in all stackings')

    print('\n5. ELECTRONIC STRUCTURE:')
    print('   - Homogeneous element (pure phosphorus) with covalent bonding')
    print('   - No inherent charge asymmetry between atoms')
    print('   - Symmetric charge distribution prevents polarization')

print('\n' + '='*70)
print('ANALYSIS COMPLETE')
print('='*70)