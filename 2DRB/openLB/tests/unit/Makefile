OLB_ROOT := ../..
include ../default.mk

TEST_OBJ_DIR := build/obj
TEST_DEP_DIR := build/dep

#TEST_GROUPS := \
#	advectionDiffusion \
#	converter \
#	functors \
#	geometry \
#	helper \
#	lattice \
#	momenta \
#	rtlbm \
#	serializable \
#	trt \
#	utilities \
#	vtk \
#	mrt \
#	postprocessor \
#	boundary \
#	dynamics \
#	particles

#TEST_CPP_FILES := $(shell find $(TEST_GROUPS) -name 'test*.cpp')

TEST_CPP_FILES := $(shell find . -name 'test*.cpp')

TEST_TARGETS := $(patsubst %.cpp, %, $(TEST_CPP_FILES))

TESTS := $(notdir $(basename $(TEST_TARGETS)))
RUNTESTS := $(addprefix run-,$(notdir $(basename $(TEST_TARGETS))))

all: $(RUNTESTS)

build: $(TESTS)

SHARED_TEST_HEADER := ./helper/unitTestHelper.h
CLEAN_TARGETS :=

CXXFLAGS += -I.

define runTest
$(1) : $(2)
	$(2)

.PHONY: $(1)
endef

define compileTest
$(TEST_OBJ_DIR)/$(1).o : $(1).cpp $(SHARED_TEST_HEADER)
	@echo $$@
	mkdir -p $$(dir $$@)
	$(CXX) $(CXXFLAGS) -DPATH='"$$(realpath $$(dir $$<))/"' -c $$< -o $$@

$(TEST_DEP_DIR)/$(1).d : $(1).cpp
	@mkdir -p $$(dir $$@)
	@sh -ec '$(CXX) -M $(CXXFLAGS) $$< \
	    | sed -e "s!$(notdir $(1)).o!$(TEST_OBJ_DIR)/$(1).o!1" > $$@;'

include $(wildcard $(TEST_DEP_DIR)/$(1).d)
endef

define targetTest
$(1) : $(2)
endef

define linkTest
$(1) : $(TEST_OBJ_DIR)/$(1).o $(TEST_DEP_DIR)/$(1).d
	$(CXX) $$< $(LDFLAGS) -o $(1)

clean-$(1):
	@rm -f $(1) &> /dev/null || true
	@rm -f $(TEST_OBJ_DIR)/$(1).o &> /dev/null || true
	@rm -f $(TEST_DEP_DIR)/$(1).d &> /dev/null || true

CLEAN_TARGETS += clean-$(1)
endef

$(foreach test,$(TEST_TARGETS), \
    $(eval $(call linkTest,$(test))))

$(foreach test,$(TEST_TARGETS), \
    $(eval $(call compileTest,$(test))))

$(foreach test,$(TEST_TARGETS), \
    $(eval $(call targetTest,$(notdir $(basename $(test))),$(test))))

$(foreach test,$(TEST_TARGETS), \
    $(eval $(call runTest,run-$(notdir $(basename $(test))),$(test))))

clean: $(CLEAN_TARGETS)
	@rm -f $(TEST_TARGETS) &> /dev/null || true

.PHONY: clean
