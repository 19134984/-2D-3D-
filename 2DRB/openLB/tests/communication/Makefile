EXAMPLE = testCommIntegrity
OLB_ROOT := ../..
override FEATURES += CASE_SILENT DIRECTIONAL_CANARY
# Expose the test-local POPULATION2 field (twoPop.h) to the mixed-mode olbcuda.cu
OLB_MIXED_MODE_INCLUDE_LOCAL_HEADERS := ON
include $(OLB_ROOT)/default.mk
