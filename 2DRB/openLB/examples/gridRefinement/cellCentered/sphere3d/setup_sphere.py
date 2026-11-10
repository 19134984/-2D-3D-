#!/usr/bin/env uv run
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy",
#   "networkx",
#   "pyvista",
#   "numba",
# ]
# ///

import sys
import os
import numpy as np

# Add the repository root to the python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..')))

from script.meshing import MesherCLI, Domain
from script.meshing.geometry import CuboidIndicator

cli = MesherCLI(default_ranks=1, default_output="sphere_setup.xml")
cli.add_case_parameter("PHYS_CHAR_LENGTH", 0.2)
cli.add_case_parameter("RESOLUTION", 16)
args = cli.parse_args()

physDeltaX = args.PHYS_CHAR_LENGTH / args.RESOLUTION
ref_block_size = 8

# Domain origin and extent
origin = np.array([0.0, 0.0, 0.0])
extent = np.array([2.5, 0.5, 0.5])

# Define Geometries / Indicators
domain_box = CuboidIndicator(
    origin=origin,
    extent=extent,
    name="domain_box"
)

domain = Domain(
    indicator=domain_box,
    base_deltaR=physDeltaX,
    root_block_size=ref_block_size
)

# Sphere bounding box refinement:
# sphere center: (0.55, 0.25, 0.25), radius: 0.1
# bounds: min=(0.45, 0.15, 0.15), max=(0.65, 0.35, 0.35)
sphere_bbox = CuboidIndicator(
    origin=[0.45, 0.15, 0.15],
    extent=[0.2, 0.2, 0.2],
    name="sphere_bbox"
)

# Trailing near field refinement (Wake region behind the sphere)
sphere_wake = CuboidIndicator(
    origin=[0.65, 0.15, 0.15],
    extent=[0.5, 0.2, 0.2],
    name="sphere_wake"
)

# Apply Refinement Rules
domain.refine(indicator=sphere_bbox, level=3, type="bounding-box", padding=6)
domain.refine(indicator=sphere_wake, level=2, type="bounding-box", padding=3)

# Generate and Export
cli.execute(domain)
