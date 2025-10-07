#!/usr/bin/env python3
"""
Create DFT-ready h-BN bilayer structures with proper cell parameters and vacuum
"""

import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import make_supercell

def create_hBN_bilayer_for_dft(stacking_type='AA', supercell_size=(1,1,1), vacuum=20.0):
    """
    Create DFT-ready h-BN bilayer structure.

    Parameters:
    - stacking_type: 'AA', 'AB', or 'BA'
    - supercell_size: tuple for supercell expansion
    - vacuum: vacuum space in Angstroms (total, split above and below)
    """

    # h-BN lattice parameters (experimental values)
    a = 2.504  # Angstroms (slightly adjusted for better DFT)
    c_spacing = 3.33  # Interlayer distance for h-BN (experimental)

    # Create hexagonal unit cell
    cell = [[a, 0, 0],
            [-a/2, a*np.sqrt(3)/2, 0],
            [0, 0, c_spacing + vacuum]]

    # Define stacking vectors (in fractional coordinates)
    tau_vectors = {
        'AA': [0.0, 0.0],
        'AB': [1/3, 2/3],  # Note: Using more standard AB stacking
        'BA': [2/3, 1/3],  # Opposite of AB
    }

    tau = tau_vectors.get(stacking_type, [0.0, 0.0])

    # Layer 1 positions (fractional coordinates)
    layer1_positions = [
        [0.0, 0.0, 0.5 - c_spacing/(2*(c_spacing + vacuum))],  # B
        [1/3, 2/3, 0.5 - c_spacing/(2*(c_spacing + vacuum))],  # N
    ]

    # Layer 2 positions (shifted by tau and raised by c_spacing)
    layer2_positions = [
        [tau[0], tau[1], 0.5 + c_spacing/(2*(c_spacing + vacuum))],  # B
        [(1/3 + tau[0]) % 1.0, (2/3 + tau[1]) % 1.0, 0.5 + c_spacing/(2*(c_spacing + vacuum))],  # N
    ]

    # Combine positions
    positions = layer1_positions + layer2_positions

    # Create atoms object
    atoms = Atoms('B2N2',
                  scaled_positions=positions,
                  cell=cell,
                  pbc=[True, True, True])

    # Create supercell if requested
    if supercell_size != (1, 1, 1):
        P = np.diag(supercell_size)
        atoms = make_supercell(atoms, P)

    return atoms


def create_all_hBN_structures():
    """Create all three h-BN stacking configurations for DFT"""

    print("Creating DFT-ready h-BN bilayer structures...")
    print("=" * 60)

    # Configuration details
    configs = {
        'AA': {
            'description': 'AA\' (bulk) stacking - Non-polar',
            'expected': 'Non-polar, high-symmetry configuration'
        },
        'AB': {
            'description': 'AB stacking - Polar (downward)',
            'expected': 'Polar with downward polarization (-P)'
        },
        'BA': {
            'description': 'BA stacking - Polar (upward)',
            'expected': 'Polar with upward polarization (+P)'
        }
    }

    # Create each configuration
    for stacking in ['AA', 'AB', 'BA']:
        print(f"\nCreating {stacking} configuration:")
        print("-" * 40)
        print(f"Description: {configs[stacking]['description']}")
        print(f"Expected: {configs[stacking]['expected']}")

        # Create structure
        atoms = create_hBN_bilayer_for_dft(
            stacking_type=stacking,
            supercell_size=(3, 3, 1),  # 3x3 supercell for better k-point sampling
            vacuum=20.0  # 20 Å vacuum
        )

        # Save in multiple formats for compatibility
        filename_cif = f'hBN_{stacking}_DFT.cif'
        filename_vasp = f'hBN_{stacking}_POSCAR'
        filename_xyz = f'hBN_{stacking}_DFT.xyz'

        write(filename_cif, atoms)
        write(filename_vasp, atoms, format='vasp', vasp5=True, direct=True)
        write(filename_xyz, atoms)

        print(f"  Cell parameters:")
        print(f"    a = {atoms.cell[0,0]:.3f} Å")
        print(f"    b = {atoms.cell[1,1]:.3f} Å")
        print(f"    c = {atoms.cell[2,2]:.3f} Å")
        print(f"  Number of atoms: {len(atoms)}")
        print(f"  Files created:")
        print(f"    - {filename_cif} (CIF format)")
        print(f"    - {filename_vasp} (VASP POSCAR format)")
        print(f"    - {filename_xyz} (XYZ format)")

        # Print recommended DFT settings
        if stacking == 'AA':
            print("\n  Recommended DFT settings:")
            print("    - Functional: PBE-D3 or optB88-vdW for van der Waals")
            print("    - K-points: 6x6x1 or higher (Gamma-centered)")
            print("    - Cutoff: 500 eV (VASP) or 80 Ry (QE)")
            print("    - Force convergence: 0.01 eV/Å")

    print("\n" + "=" * 60)
    print("DFT-ready structures created successfully!")
    print("\nNotes for DFT calculations:")
    print("1. The structures have 20 Å vacuum to prevent periodic interactions")
    print("2. 3x3 supercells are used for better k-point sampling")
    print("3. Interlayer distance set to 3.33 Å (experimental value)")
    print("4. You may need to relax the structures to find optimal geometry")
    print("5. For polarization calculations, use Berry phase methods")
    print("=" * 60)


if __name__ == "__main__":
    create_all_hBN_structures()