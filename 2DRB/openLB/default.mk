ifndef IN_OLB_NIX_SHELL
# Include config.mk environment (optional)
-include $(OLB_ROOT)/config.mk
endif

# Select mixed compilation mode if separate CUDA compiler is given
ifdef SYCL_CXX
include $(OLB_ROOT)/default.sycl.mixed.mk
else
ifdef CUDA_CXX
include $(OLB_ROOT)/default.mixed.mk
else
include $(OLB_ROOT)/default.single.mk
endif
endif
