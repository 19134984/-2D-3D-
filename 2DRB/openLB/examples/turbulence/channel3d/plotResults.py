# Copyright (C) 2025 Fedor Bukreev
# E-mail contact: info@openlb.net
# The most recent release of OpenLB can be downloaded at
# <http://www.openlb.net/>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# --- Configuration ---
computed_file = 'tmp/gnuplotData/data/flowField.csv'
reference_file = 'referenceData.csv'
output_dir = 'comparativePlots'

# Define the common x-axis columns
x_computed_col = 'y+'
x_reference_col = 'y+'

# Define the column mappings for comparison
# Format: (computed_col, reference_col, plot_title, y_label, filename_suffix)
plot_mappings = [
    ('uAv+', 'u+', 'Mean Velocity Profile', '$u^+$', 'u_plus'),
]
# ---------------------

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")
else:
    print(f"Output directory already exists: {output_dir}")

# Load Computed Data
try:
    computed_df = pd.read_csv(computed_file, delimiter=';')
    print(f"\n--- Successfully loaded computed data: {computed_file} ---")
    print("Head:")
    print(computed_df.head())
    print("\nInfo:")
    computed_df.info()
except FileNotFoundError:
    print(f"\nError: Computed data file not found at {computed_file}")
    print(f"Please make sure '{computed_file}' is uploaded.")
    sys.exit(1) # Exit script if a file is missing
except Exception as e:
    print(f"\nError loading {computed_file}: {e}")
    sys.exit(1)

# Load Reference Data
try:
    reference_df = pd.read_csv(reference_file, delimiter=',')
    print(f"\n--- Successfully loaded reference data: {reference_file} ---")
    print("Head:")
    print(reference_df.head())
    print("\nInfo:")
    reference_df.info()
except FileNotFoundError:
    print(f"\nError: Reference data file not found at {reference_file}")
    print(f"Please make sure '{reference_file}' is uploaded.")
    sys.exit(1)
except Exception as e:
    print(f"\nError loading {reference_file}: {e}")
    sys.exit(1)

# Generate Plots
print("\n--- Generating plots ---")
generated_files = []

for comp_col, ref_col, title, y_label, suffix in plot_mappings:
    # Check if all required columns exist
    if comp_col not in computed_df.columns:
        print(f"Warning: Column '{comp_col}' not found in {computed_file}. Skipping plot '{title}'.")
        continue
    if x_computed_col not in computed_df.columns:
        print(f"Warning: Column '{x_computed_col}' not found in {computed_file}. Skipping plot '{title}'.")
        continue
    if ref_col not in reference_df.columns:
        print(f"Warning: Column '{ref_col}' not in {reference_file}. Skipping plot '{title}'.")
        continue
    if x_reference_col not in reference_df.columns:
        print(f"Warning: Column '{x_reference_col}' not in {reference_file}. Skipping plot '{title}'.")
        continue

    try:
        plt.figure(figsize=(10, 6))

        plt.plot(computed_df[x_computed_col], computed_df[comp_col],
                 marker='o', linestyle='-', color='blue',
                 label='Computed', markersize=4, markerfacecolor='none')
        plt.plot(reference_df[x_reference_col], reference_df[ref_col],
                 linestyle='--', color='red',
                 label='Reference')

        plt.xscale('log')  # Use logarithmic scale for y+ axis
        plt.title(f"{title} ({y_label}) vs. $y^+$")
        plt.xlabel('$y^+$')
        plt.ylabel(y_label)
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.6)
        plt.tight_layout()

        filename = os.path.join(output_dir, f'comparison_{suffix}.png')
        plt.savefig(filename)
        plt.close()  # Close the figure to free memory
        generated_files.append(filename)
        print(f"Successfully generated: {filename}")

    except Exception as e:
        print(f"Error generating plot for '{title}': {e}")
        plt.close() # Ensure figure is closed even if an error occurs

print(f"\n--- Plot generation complete ---")
print(f"Total plots generated: {len(generated_files)}")
print("Generated files:")
for f in generated_files:
    print(f)
