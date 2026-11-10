# Build config for Intel oneAPI C++ (icpx) continuous-integration coverage.
#
# Plain CPU compilation with the Intel LLVM-based compiler, no SYCL offload.
# Runs inside the stock public oneAPI hpckit image (intel/oneapi-hpckit), whose
# system gcc supplies the standard library, so no --gcc-toolchain is needed as
# on HoreKa.

CXX             := icpx
CC              := icx

CXXFLAGS        := -O3 -Wall
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD

USE_EMBEDDED_DEPENDENCIES := ON
