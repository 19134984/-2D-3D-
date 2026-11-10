# OpenLB build configuration
#
# This file sets up the necessary build flags for compiling OpenLB with
# the GNU C++ compiler and sequential execution. For more complex setups
# edit this file or consult the example configs provided in `config/`.
#
# Basic usage:
#  - Edit variables to fit desired configuration
#  - Run `make clean; make` to clean up any previous artifacts and compile the dependencies
#  - Switch to example directory, e.g. `examples/laminar/poiseuille2d`
#  - Run `make`
#  - Start the simulation using `./poiseuille2d`

# Compiler to use for C++ files, change to `mpic++` when using OpenMPI and GCC
CXX             := g++
# Compiler to use for C files (used for emebedded dependencies)
CC              := gcc

# Suggested optimized build flags for GCC, consult `config/` for further examples
# for clang use -march=native
CXXFLAGS        := -O3 -Wall -march=native -mtune=native
# Uncomment to add debug symbols and enable runtime asserts
#CXXFLAGS        += -g -DOLB_DEBUG

# OpenLB requires support for C++20
# works in:
#  * gcc 12 or later     (https://gcc.gnu.org/projects/cxx-status.html#cxx20)
#  * clang 19 or later   (https://clang.llvm.org/cxx_status.html#cxx20)
#  * icx 2024.0 or later (https://www.intel.com/content/www/us/en/developer/articles/technical/c20-features-supported-by-intel-cpp-compiler.html)
CXXFLAGS        += -std=c++20

# optional linker flags
LDFLAGS         :=

# Parallelization mode, must be one of: OFF, MPI, OMP, HYBRID
# Note that for MPI and HYBRID the compiler also needs to be adapted.
# See e.g. `config/cpu_gcc_openmpi.mk`
PARALLEL_MODE   := OFF

# optional MPI and OpenMP flags
MPIFLAGS        :=
OMPFLAGS        := -fopenmp

# Options: CPU_SISD, CPU_SIMD, GPU_CUDA
# Both CPU_SIMD and GPU_CUDA require system-specific adjustment of compiler flags.
# See e.g. `config/cpu_simd_intel_mpi.mk` or `config/gpu_only.mk` for examples.
# CPU_SISD must always be present.
PLATFORMS       := CPU_SISD

# Fundamental arithmetic data type
# Common options are float or double
FLOATING_POINT_TYPE := double

# Optional features, passed to the compiler as `-DFEATURE_<name>` declarations.
# Enable by listing space-separated names, e.g. `FEATURES := VTK ADIOS2`.
# Recognized features include:
#   Dependencies / I/O:  VTK, ADIOS2, VDB, PROJ, OPENBLAS, NLOPT
#   Diagnostics:         PROFILING, INSPECT_DYNAMICS, INSPECT_POST_PROCESSORS, HIDE_LOGO
#   Code generation:     DISABLE_CSE (recompute collisions instead of using CSE)
#   GPU / portability:
#     <PLATFORM>_STORED_MOMENTS
#                           Store rho/u (and pi/j) as fields during collision instead
#                           of recomputing them via virtual dispatch. The demand
#                           resolver then refreshes boundary neighbours for post-stream
#                           freshness. Whole-lattice and coupling reads see post-collide
#                           values, so keep the setting uniform across platforms for
#                           cross-platform reproducibility. Useful where couplings read
#                           moments over the whole bulk. See src/core/platform/platform.h.
#                             | Platform | Feature flag            | Default |
#                             |----------+-------------------------+---------|
#                             | GPU_SYCL | always on (structural)  | on      |
#                             | GPU_HIP  | HIP_STORED_MOMENTS      | off     |
#                             | GPU_CUDA | CUDA_STORED_MOMENTS     | off     |
#                             | CPU_SIMD | CPU_SIMD_STORED_MOMENTS | off     |
#                             | CPU_SISD | CPU_SISD_STORED_MOMENTS | off     |
#     FSI_ATOMIC_REDUCTION  Use atomicAdd scatter for the GPU FSI element
#                           reduction instead of the default sort + reduce_by_key.
#                           See src/fsi/platform/gpu/{cuda,hip}/integral.hh.
#     GPU_CYCLIC_BRANCHING  Force the branching two-pointer cyclic population store.
#     GPU_CYCLIC_VMEM       Force the virtual-memory ring cyclic population store.
#                           (Default: vmem on CUDA, branching on HIP.)
#                           See src/core/platform/gpu/cyclicAccessor.h.
FEATURES        :=

# The HIP target architecture is set via the HIP_ARCH make variable (e.g. gfx90a); the
# nix `hip`/`hip-cuda` dev shells take it from the OLB_HIP_ARCH environment variable.

# Set to OFF if libz and tinyxml2 are provided by the system (optional)
USE_EMBEDDED_DEPENDENCIES := ON
