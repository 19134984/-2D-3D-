include $(OLB_ROOT)/rules.mk

export LC_MESSAGES=C

all: dependencies core $(EXAMPLE)

INCLUDE_DIRS := $(OLB_ROOT)/src

ifneq ($(USE_EMBEDDED_DEPENDENCIES), OFF)
INCLUDE_DIRS += \
	$(OLB_ROOT)/external/zlib \
	$(OLB_ROOT)/external/tinyxml2

LDFLAGS := -L$(OLB_ROOT)/external/lib $(LDFLAGS)

dependencies:
	$(MAKE) -C $(OLB_ROOT)/external

clean-dependencies:
	$(MAKE) -C $(OLB_ROOT)/external clean
else
.PHONY: dependencies clean-dependencies
endif

LDFLAGS += -L$(OLB_ROOT)/build/lib

core:
	$(MAKE) -C $(OLB_ROOT) core

CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)
DEP_FILES := $(CPP_FILES:.cpp=.d)

INCLUDE_FLAGS := -I$(subst $(EMPTY) $(EMPTY), -I,$(INCLUDE_DIRS))

%.d: $(CPP_FILES)
	@$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDE_FLAGS) $< -MM -MT $(@:.d=.o) >$@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

build/missing.txt: $(OBJ_FILES)
	mkdir -p build
	$(CXX) $^ $(LDFLAGS) -lolbcore 2>&1 | c++filt --no-recurse-limit | grep -oP "undefined reference to \`\K[^']+\)" | sort | uniq > $@

EXPLICIT_METHOD_INSTANTIATION := \
	FieldTypeRegistry \
	StatisticsPostProcessor \
	FieldStatisticO \
	IntegratePorousElementFieldsO \
	BlockLatticeFieldReduction \
	checkPlatform \
	getFusedCollision

EXPLICIT_CLASS_INSTANTIATION := \
	Column \
	CyclicColumn \
	ConcreteBlockO \
	ConcreteBlockCollisionO \
	ConcreteBlockPointCouplingO \
	ConcreteBlockRefinementO \
	ConcreteCommunicatable \
	ConcreteBlockCommunicator \
	HeterogeneousCopyTask \
	ConcreteParametersD \
	BlockLatticeFieldReduction \
	unique_ptr

EXPLICIT_CLASS_INSTANTIATION_LINKER := \
	Column \
	CyclicColumn \
	ConcreteBlockO \
	ConcreteBlockCollisionO \
	ConcreteBlockCouplingO \
	ConcreteBlockPointCouplingO \
	ConcreteBlockRefinementO \
	ConcreteCommunicatable \
	ConcreteBlockCommunicator \
	HeterogeneousCopyTask \
	ConcreteParametersD \
	unique_ptr \
	StatisticsPostProcessor \
	FieldStatisticO \
	IntegratePorousElementFieldsO \
	BlockLatticeFieldReduction \
	Event

CPP_HASH := \#
OLB_MIXED_MODE_INCLUDE := ${CPP_HASH}include <olb.h>\n
ifdef OLB_MIXED_MODE_INCLUDE_CPP
    OLB_MIXED_MODE_INCLUDE := $(CPP_FILES:%=\n${CPP_HASH}include "../%"\n)
endif

ifeq ($(OLB_MIXED_MODE_INCLUDE_LOCAL_HEADERS), ON)
    LOCAL_HEADER_FILES := $(wildcard *.h)
    OLB_MIXED_MODE_INCLUDE += $(LOCAL_HEADER_FILES:%=\n${CPP_HASH}include "../%"\n)
endif

build/olbsycl.sycl: build/missing.txt
	printf '$(OLB_MIXED_MODE_INCLUDE)' > $@
# Transform missing symbols into explicit template instantiations by:
# - filtering for a set of known and automatically instantiable methods
# - excluding destructors
# - dropping resulting empty lines
# - adding the explicit instantiation prefix (all supported methods are void, luckily)
	cat build/missing.txt \
	| grep '$(subst $() $(),\|,$(EXPLICIT_METHOD_INSTANTIATION))' \
	| grep -wv '.*~.*\|FieldTypeRegistry()' \
	| grep -v 'FieldTypeRegistry<.*(olb::Platform)[0-9]>::[[:alnum:]_]\+(' \
	| sed -e 's/.*/template void &;/' -e 's/void void/void/' -e 's/void std::function/std::function/' >> $@
# | xargs -0 -n1 | grep . \
# - filtering for a set of known and automatically instantiable classes
# - dropping method cruft and wrapping into explicit class instantiation
# - removing duplicates
	cat build/missing.txt \
	| grep '.*\($(subst $() $(),\|,$(EXPLICIT_CLASS_INSTANTIATION))\)<' \
	| sed -e 's/>::~\?[[:alnum:]]\+\[\?\]\?(.*)/>/' -e 's/.*/template class &;/' -e 's/class void/class/' \
	| sort | uniq >> $@
# Handle ConcreteBlockCouplingO
	cat build/missing.txt \
	| grep '.*\($(subst $() $(),\|,ConcreteBlockCouplingO)\)<' \
	| grep -v 'ConcreteBlockCouplingO.*ConcreteBlockCouplingO' \
	| sed -e 's/>::~\?[[:alnum:]]\+\[\?\]\?(.*)/>/' -e 's/.*/template class &;/' -e 's/class void/class/' \
	| sort | uniq >> $@
	cat build/missing.txt \
	| sed -n 's/\(olb::ConcreteBlockCouplingO.*ConcreteBlockCouplingO.*\)/\1/p' \
	| sed -e 's/>::ConcreteBlockCouplingO<.*>(/>::ConcreteBlockCouplingO(/' -e 's/.*/template &;/' \
	| sort | uniq >> $@
# Handle FieldTypeRegistry
	cat build/missing.txt \
  | sed -n 's/\(olb::FieldTypeRegistry<.*(olb::Platform)[0-9]>\)::.*/\1/p' \
	| sed -e 's/.*/template class &;/' -e 's/class void/class/' \
	| sort | uniq >> $@

SYCL_CPP_FILES := build/olbsycl.sycl
SYCL_OBJ_FILES := $(SYCL_CPP_FILES:.sycl=.o)

ifneq ($(strip $(GIT_VERSION)),)
	SYCL_CXXFLAGS += -DOLB_VERSION=\"$(GIT_VERSION)\"
else
	SYCL_CXXFLAGS += -DOLB_VERSION=\"$(OLB_RELEASE)\"
endif

SYCL_CXXFLAGS += $(filter -DFEATURE_%,$(CXXFLAGS))

%.o: %.sycl
	$(SYCL_CXX) -x c++ $(SYCL_CXXFLAGS) $(INCLUDE_FLAGS) $(PARALLEL_FLAGS) -DPLATFORM_CPU_SISD -DPLATFORM_GPU_SYCL -DOLB_PRECOMPILING_SYCL -fPIC -c -o $@ $<

build/olbsycl.version: $(SYCL_OBJ_FILES)
	echo 'LIBOLBSYCL { global: ' > $@
# Declare exposed explicitly instantiated symbols to prevent duplicate definitions by:
# - filtering for the set of automatically instantiated classes
# - excluding CPU_SISD symbols (we only instantiate GPU_SYCL-related symbols)
# - dropping the shared library location information
# - postfixing by semicolons
	nm $(SYCL_OBJ_FILES) \
	| grep '$(subst $() $(),\|,$(EXPLICIT_CLASS_INSTANTIATION_LINKER))\|sycl.*device\|checkPlatform\|getFusedCollision\|FieldTypeRegistry' \
	| grep -wv '.*sisd.*' \
	| cut -c 20- \
	| sed 's/$$/;/' >> $@
	echo 'local: *; };' >> $@

libolbsycl.so: $(SYCL_OBJ_FILES) build/olbsycl.version
	$(SYCL_CXX) $(SYCL_CXXFLAGS) -Xlinker --version-script=build/olbsycl.version -shared $(SYCL_OBJ_FILES) -o $@

$(EXAMPLE): $(OBJ_FILES) libolbsycl.so
	$(CXX) $(OBJ_FILES) -o $@ $(LDFLAGS) -L . -lolbsycl -lolbcore $(SYCL_LDFLAGS)

$(EXAMPLE)-no-sycl-recompile: $(OBJ_FILES)
	$(CXX) $^ -o $(EXAMPLE) $(LDFLAGS) -L . -lolbsycl -lolbcore $(SYCL_LDFLAGS)

.PHONY: onlysample
onlysample: $(EXAMPLE)

.PHONY: no-sycl-recompile
no-sycl-recompile: $(EXAMPLE)-no-sycl-recompile

.PHONY: clean-tmp
clean-tmp:
	rm -f tmp/*.* tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/imageData/data/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

.PHONY: clean-core
clean-core:
	make -C $(OLB_ROOT) clean-core

.PHONY: clean
clean: clean-tmp
	rm -f $(OBJ_FILES) $(DEP_FILES) $(EXAMPLE) $(SYCL_OBJ_FILES) $(SYCL_CPP_FILES) libolbsycl.so

-include $(DEP_FILES)
