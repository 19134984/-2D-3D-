CXX             := icpx
CC              := icx

CXXFLAGS := -O3
CXXFLAGS += -std=c++20

CXXFLAGS += -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc -O0"
LDFLAGS  += -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc"

# intermediate output and full log for debugging
LDFLAGS += -save-temps -v

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD GPU_SYCL

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
