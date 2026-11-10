# Build config for GPU_SYCL continuous-integration coverage.
#
# Targets the SYCL backend with the Intel oneAPI DPC++ compiler (icpx) without a
# fixed device target, so the JIT-compiled SPIR-V kernels run on whichever SYCL
# device is selected at runtime. On a CI runner without a GPU this is the OpenCL
# CPU device (select via ONEAPI_DEVICE_SELECTOR=opencl:cpu), which makes the
# resulting binaries runnable for smoke testing.
#
# Usage (inside an Intel oneAPI hpckit container):
#  - Copy this file to OpenLB root as `config.mk`
#  - Run `make`

CXX             := icpx
CC              := icx

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD GPU_SYCL

SYCL_CXX        := icpx
SYCL_CXXFLAGS   += -O3 -std=c++20 -fsycl -fsycl-device-code-split=per_kernel
SYCL_LDFLAGS    += -fsycl -fsycl-device-code-split=per_kernel

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
