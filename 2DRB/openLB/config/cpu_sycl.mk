CXX             := mpicxx
CC              := icx

CXXFLAGS := -O3 -g -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_SYCL

SYCL_CXX := icpx
SYCL_CXXFLAGS += -O3 -std=c++20 -fsycl -fsycl-device-code-split=per_kernel
SYCL_LDFLAGS  += -fsycl -fsycl-device-code-split=per_kernel

FEATURES := VTK
VTK_VERSION := 9.1

CXXFLAGS += -I/usr/include/vtk-9.1

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
