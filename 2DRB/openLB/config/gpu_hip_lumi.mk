CXX             := CC

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++20 -x hip 
CXXFLAGS        += -I$(ROCM_PATH)/include -I$(ROCM_PATH)/include/hip
PARALLEL_MODE   := MPI

HIP_PLATFORM	:= amd  
PLATFORMS       := CPU_SISD GPU_HIP

HIP_ARCH	:= gfx90a #mi250

FEATURES        := ADIOS2

FLOATING_POINT_TYPE := float

