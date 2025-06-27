#!/usr/bin/env python3
"""
Runner script for SYMPOL2D
"""

import sys
import os

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sympol2d import cli

if __name__ == '__main__':
    sys.exit(cli.main())