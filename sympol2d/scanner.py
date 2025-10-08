"""
Stacking configuration scanner for SYMPOL2D
"""

import numpy as np
from .symmetry import LAYER_GROUP_OPERATIONS, survives, classify_z_allowed, Z_SIGN_FLIP_EXPECTED, mod1


def grid_points(n):
    """Generate n×n grid of fractional coordinates in [0,1)"""
    xs = np.linspace(0, 1, n, endpoint=False)
    return np.array([(x, y) for x in xs for y in xs], float)


def survivors_at_tau(layer_group, tau, atol=1e-8):
    """
    Get list of symmetry operations that survive at stacking vector tau.

    Args:
        layer_group: Layer group symbol (e.g., 'p-6m2', 'pman')
        tau: Stacking vector [tau_x, tau_y] in fractional coordinates
        atol: Tolerance for integer check

    Returns:
        List of operation names that survive
    """
    ops = LAYER_GROUP_OPERATIONS.get(layer_group)
    if not ops:
        raise ValueError(f"Unsupported layer group: {layer_group}")
    return [R for R in ops if survives(R, np.array(tau, float), atol=atol)]


def scan_for_z(layer_group, ngrid=60, atol=1e-8):
    """
    Scan tau-space for z-polarization allowed stackings.

    Args:
        layer_group: Layer group symbol
        ngrid: Grid resolution (ngrid × ngrid points)
        atol: Tolerance for symmetry test

    Returns:
        List of candidate dictionaries with keys:
            - tau: stacking vector
            - survivors: list of surviving operations
            - tag: "z-only" or "z-allowed"
    """
    pts = grid_points(ngrid)
    keep = []
    for t in pts:
        surv = survivors_at_tau(layer_group, t, atol=atol)
        z_ok, tag = classify_z_allowed(surv)
        if z_ok:
            keep.append({"tau": t, "survivors": surv, "tag": tag})
    return keep


def distance_to_mirror_lines(t):
    """
    Calculate distance to nearest mirror lines.
    Mirrors at tau_x in {0, 0.5}, tau_y in {0, 0.5}

    Args:
        t: Stacking vector [tau_x, tau_y]

    Returns:
        Sum of distances to nearest mirror lines
    """
    dx = min(abs(t[0]-0.0), abs(t[0]-0.5), abs(t[0]-1.0))
    dy = min(abs(t[1]-0.0), abs(t[1]-0.5), abs(t[1]-1.0))
    return dx + dy


def rational_score(t):
    """
    Score for preferring simple rational fractions.

    Args:
        t: Stacking vector [tau_x, tau_y]

    Returns:
        Score (lower is better) based on distance to simple fractions
    """
    fracs = np.array([1/2, 1/3, 2/3, 1/4, 3/4, 1/6, 5/6])
    return np.min(abs(t[0]-fracs)) + np.min(abs(t[1]-fracs))


def pick_best_pair(layer_group, candidates, prefer_strict=True):
    """
    Select exactly ONE best (tau, 1-tau) pair from candidates.

    Selection criteria (in order):
        1) Prefer 'z-only' tags (mirrors broken, C2 survives)
        2) Maximize distance from mirror lines
        3) Prefer simple rational fractions

    Args:
        layer_group: Layer group symbol
        candidates: List of candidate dictionaries from scan_for_z
        prefer_strict: If True, prefer "z-only" over "z-allowed"

    Returns:
        (tau, tau_ba): Best AB/BA pair, or (None, None) if no candidates
    """
    cands = candidates
    if prefer_strict:
        strict = [c for c in cands if c["tag"] == "z-only"]
        if strict:
            cands = strict

    def key(c):
        t = c["tau"]
        # Prioritize simple rational fractions first, then maximize distance from mirrors
        return (rational_score(t), -distance_to_mirror_lines(t))

    cands = sorted(cands, key=key)
    if not cands:
        return None, None
    tau = cands[0]["tau"]
    return tau, mod1(1.0 - tau)


def z_sign_flip_expected(layer_group):
    """
    Check if AB↔BA sliding is expected to flip Pz sign for this layer group.

    Args:
        layer_group: Layer group symbol

    Returns:
        bool: True if Pz flip expected (hexagonal families),
              False if not (rectangular/square families)
    """
    return bool(Z_SIGN_FLIP_EXPECTED.get(layer_group, False))
