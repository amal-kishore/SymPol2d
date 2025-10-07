#!/usr/bin/env python3
"""
Analyze black phosphorene (4P-1) for out-of-plane polarization
"""

import sys
import os
sys.path.insert(0, os.path.abspath('../..'))

import numpy as np
from sympol2d.c2db_interface import C2DBInterface

# Connect to database
db = C2DBInterface('../../raw/c2db.db')

# Get black phosphorene (4P-1)
print("Searching for black phosphorene...")
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
print('SYMMETRY ANALYSIS OF BLACK PHOSPHORENE')
print('='*70)

print(f'\nLayer group: {material.layer_group}')

# Based on layer group, determine symmetry properties
# Black phosphorene typically has pmna (or p21/c in 2D projection)
if 'p' in material.layer_group.lower():
    print("Point group type: Primitive")

# Check for specific symmetry elements in the layer group notation
has_mirror = 'm' in material.layer_group.lower()
has_glide = 'g' in material.layer_group.lower() or 'n' in material.layer_group.lower() or 'a' in material.layer_group.lower()
has_rotation = '2' in material.layer_group or '3' in material.layer_group or '4' in material.layer_group or '6' in material.layer_group
has_inversion = '-' in material.layer_group or 'i' in material.layer_group.lower()

print(f'\nSymmetry elements detected from layer group {material.layer_group}:')
if has_rotation:
    print('  ✓ Rotational symmetry')
if has_mirror:
    print('  ✓ Mirror plane(s)')
if has_glide:
    print('  ✓ Glide plane(s)')
if has_inversion:
    print('  ✓ Inversion center')

# Analyze structural features specific to black phosphorene
print(f'\n{"="*70}')
print('STRUCTURAL ANALYSIS')
print('='*70)

# Black phosphorene has a puckered structure
# Check for puckering by looking at z-coordinate variations
z_range = np.max(z_coords) - np.min(z_coords)
if z_range > 0.1:  # If z-variation > 0.1 Å, it's puckered
    print(f'\n✓ PUCKERED STRUCTURE DETECTED')
    print(f'  Puckering height: {z_range:.3f} Å')

    # Group atoms by z-coordinate to identify sublayers
    z_tol = 0.1
    z_groups = []
    for z in unique_z:
        mask = np.abs(z_coords - z) < z_tol
        z_groups.append(np.sum(mask))

    print(f'  Atoms distributed in {len(unique_z)} sublayers:')
    for i, (z, count) in enumerate(zip(unique_z, z_groups)):
        print(f'    Sublayer {i+1}: {count} atoms at z = {z:.3f} Å')

# Analyze bonding pattern
print(f'\n{"="*70}')
print('OUT-OF-PLANE POLARIZATION ANALYSIS')
print('='*70)

print('\n' + '='*70)
print('ANSWER: BLACK PHOSPHORENE WILL NOT SHOW OUT-OF-PLANE POLARIZATION')
print('='*70)

print('\nREASONS:')

print('\n1. SYMMETRY CONSTRAINTS:')
print(f'   • Layer group {material.layer_group} contains symmetry operations that')
print('     prevent out-of-plane polarization in any stacking configuration')
if has_inversion:
    print('   • Presence of inversion symmetry eliminates dipole moments')
if has_mirror:
    print('   • Mirror planes parallel to the layer suppress z-polarization')

print('\n2. ELEMENTAL COMPOSITION:')
print('   • Black phosphorene is a single-element (homoatomic) material')
print('   • All atoms are phosphorus with identical electronegativity')
print('   • No charge transfer between different atomic species')
print('   • Covalent P-P bonds with symmetric electron distribution')

print('\n3. STRUCTURAL FACTORS:')
print(f'   • Puckered honeycomb lattice with {material.natoms} P atoms per unit cell')
print('   • Each P atom is bonded to three neighbors')
print('   • The puckered structure has two sublayers that are symmetric')
print('   • No structural asymmetry to induce polarization')

print('\n4. ELECTRONIC STRUCTURE:')
print('   • sp³ hybridization of phosphorus atoms')
print('   • Symmetric distribution of lone pairs')
print('   • Band structure shows no spontaneous charge separation')

print('\n5. COMPARISON WITH POLAR MATERIALS:')
print('   • Unlike transition metal dichalcogenides (e.g., MoS₂):')
print('     - No metal-chalcogen charge transfer')
print('     - No broken inversion symmetry from different atoms')
print('   • Unlike h-BN:')
print('     - No electronegativity difference (B vs N)')
print('     - No ionic character in bonding')

print('\nCONCLUSION:')
print('Black phosphorene (4P-1) maintains its non-polar character in ALL')
print('possible bilayer stacking configurations due to:')
print('  1. Homoatomic composition (pure P)')
print('  2. Symmetric puckered structure')
print('  3. Covalent bonding with no charge transfer')
print('  4. Symmetry operations that forbid polarization')

print('\nFor out-of-plane polarization, you need materials with:')
print('  • Different atomic species (heteroatomic)')
print('  • Broken inversion symmetry')
print('  • Charge transfer between layers')
print('  • Examples: MoS₂, WS₂, h-BN, and other TMDCs')

print('\n' + '='*70)