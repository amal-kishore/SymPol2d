"""
Bilayer structure builder for SYMPOL2D
"""

import numpy as np


def wrap01(x):
    """Wrap fractional coordinates to [0,1)"""
    return x - np.floor(x)


def make_bilayer_poscar(formula, lattice, frac_coords_mono, tau, dz_frac,
                        comment="bilayer", top_first=False):
    """
    Build rigid bilayer POSCAR from monolayer structure.

    Args:
        formula: Chemical formula string (e.g., "MoS2", "BN")
        lattice: 3×3 lattice vectors in Angstroms
        frac_coords_mono: N×3 fractional coordinates of monolayer atoms
        tau: [tau_x, tau_y] in-plane shift for top layer (fractional)
        dz_frac: Fractional vertical offset between layers
        comment: Comment line for POSCAR
        top_first: If True, list top-layer atoms first (for visualization)

    Returns:
        str: Complete POSCAR content
    """
    bottom = np.array(frac_coords_mono, float)
    top = np.array(frac_coords_mono, float)

    # Apply in-plane shift to top layer
    top[:, 0] = wrap01(top[:, 0] + tau[0])
    top[:, 1] = wrap01(top[:, 1] + tau[1])

    # Apply vertical shift to top layer
    top[:, 2] = wrap01(top[:, 2] + dz_frac)

    # Stack layers
    coords = np.vstack([top, bottom]) if top_first else np.vstack([bottom, top])

    # Build POSCAR content
    lines = []
    lines.append(f"{comment}")
    lines.append("1.0")
    for v in lattice:
        lines.append(f"  {v[0]:18.12f}  {v[1]:18.12f}  {v[2]:18.12f}")
    lines.append(formula)
    lines.append(str(coords.shape[0]))
    lines.append("Direct")
    for r in coords:
        lines.append(f"  {r[0]:.12f}  {r[1]:.12f}  {r[2]:.12f}")

    return "\n".join(lines) + "\n"
