"""
SYMPOL2D: SYMmetry-based prediction of POLarity in 2D bilayers

A tool for predicting polar and non-polar stacking configurations in 
2D bilayer systems using symmetry analysis without DFT calculations.
"""

__version__ = "0.1.0"
__author__ = "SYMPOL2D Development Team"

from . import cli
from . import symmetry
from . import scanner
from . import c2db_interface
from . import utils