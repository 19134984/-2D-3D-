###########################################################################
## OpenLB Makefile

DEFAULT_TARGETS :=
CLEAN_TARGETS :=

# Call DEFAULT_TARGETS for plain `make`
all: default

###########################################################################
## Default build configuration

# Use environment if in an OpenLB Nix shell, otherwise use config.mk
ifndef IN_OLB_NIX_SHELL
ifeq ($(wildcard config.mk),)
config.mk: config.git.mk
	@echo Copying default config
	cp -n config.git.mk config.mk

DEFAULT_TARGETS += config.mk
endif

# Include config.mk environment (optional)
-include config.mk
endif

include rules.mk

###########################################################################
## Embedded dependencies (optional)

ifneq ($(USE_EMBEDDED_DEPENDENCIES), OFF)
dependencies:
	$(MAKE) CXX='$(CXX)' CC='$(CC)' -C external

clean-dependencies:
	$(MAKE) CXX='$(CXX)' CC='$(CC)' -C external clean

CLEAN_TARGETS += clean-dependencies
else
.PHONY: dependencies
endif

DEFAULT_TARGETS += dependencies

###########################################################################
## Core library

CORE_CPP_FILES := \
	src/communication/mpiManager.cpp \
	src/communication/ompManager.cpp \
	src/core/olbInit.cpp \
	src/core/expr.cpp \
	src/io/ostreamManager.cpp \
	src/utilities/log/log.cpp

CORE_OBJ_FILES := $(CORE_CPP_FILES:.cpp=.o)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -fPIC -Isrc/ -c $< -o $@

build/lib:
	mkdir -p build/lib

build/lib/libolbcore.a: build/lib $(CORE_OBJ_FILES)
	$(AR) rc $@ $(CORE_OBJ_FILES)

core: build/lib/libolbcore.a

DEFAULT_TARGETS += core

clean-core:
	rm -f $(CORE_OBJ_FILES)
	rm -f build/lib/libolbcore.a

CLEAN_TARGETS += clean-core

###########################################################################
## Examples

EXAMPLES := $(dir $(shell find examples -name 'Makefile' -not -path 'examples/web/*'))

CLEAN_SAMPLES_TARGETS :=

.PHONY: $(EXAMPLES)
$(EXAMPLES): dependencies core
	$(MAKE) -C $@ onlysample
samples: $(EXAMPLES)

.PHONY: run-samples
run-samples:
	@for dir in $(EXAMPLES); do \
		$(MAKE) -C $$dir run; \
	done

define clean_directory
clean-$(1):
	+$(MAKE) -C $(1) clean

CLEAN_SAMPLES_TARGETS += clean-$(1)
endef

$(foreach sample,$(EXAMPLES), \
  $(eval $(call clean_directory,$(sample))))

clean-samples: $(CLEAN_SAMPLES_TARGETS)

###########################################################################
## Code generation (CSE)

CSE_GENERATORS := $(shell find src/ -type f -name '*.cse.h.template')
CSE_GENERATEES := $(patsubst %.cse.h.template, %.cse.h, $(CSE_GENERATORS))

%.cse.h: %.cse.h.template
	$(eval $@_GUARD := $(shell grep -Po "#ifndef\ \K[A-Z_]*_CSE_H" $<))
	python script/codegen/legacy/cse.py $< $($@_GUARD) $@

cse: $(CSE_GENERATEES)

###########################################################################
## Unit tests

check: dependencies core
	$(MAKE) -C tests/unit

clean-check:
	$(MAKE) -C tests/unit clean

run-%:
	$(MAKE) -C tests/unit $@

CLEAN_TARGETS += clean-check

###########################################################################
## Benchmark tests

benchmark: dependencies core
	$(MAKE) -C tests/benchmarks

test%:
	$(MAKE) -C tests/benchmarks run-$@

clean-benchmark:
	$(MAKE) -C tests/benchmarks clean

CLEAN_TARGETS += clean-benchmark

###########################################################################
## User guide documentation

userguide:
	$(MAKE) -C doc

###########################################################################
## Doxygen documentation

# The doxygen-awesome stylesheet is optional and supplied by the nix doxygen
# environment, so a plain checkout still builds the docs without it. The version
# always comes from VERSION rather than being repeated in the Doxyfile.
DOXYGEN_AWESOME_CSS ?=

doxygen:
	@mkdir -p build/doxygen
	@{ cat Doxyfile; \
	   echo "PROJECT_NUMBER = $(OLB_RELEASE)"; \
	   if [ -f doc/_doxygen_header.html ]; then \
	     echo "HTML_HEADER = doc/_doxygen_header.html"; \
	   fi; \
	   if [ -n "$(DOXYGEN_AWESOME_CSS)" ]; then \
	     echo "HTML_COLORSTYLE = LIGHT"; \
	     echo "HTML_EXTRA_STYLESHEET = $(DOXYGEN_AWESOME_CSS)/doxygen-awesome.css $(DOXYGEN_AWESOME_CSS)/doxygen-awesome-sidebar-only.css doc/doxygen-quarto.css"; \
	   fi; \
	 } | doxygen -

# The whole documentation set: the guide, then the API reference carrying the
# guide's navbar, which only exists once the guide has been rendered.
docs:
	$(MAKE) userguide
	python3 doc/doxygen_header.py
	$(MAKE) doxygen

###########################################################################
## Install git hooks

ifeq ($(shell test -d .git -a -d hooks; echo $$?),0)
.git/hooks:
	@echo Installing git hooks
	ln -rs hooks .git/hooks

DEFAULT_TARGETS += .git/hooks
endif

###########################################################################
## Cleaning

clean: $(CLEAN_TARGETS)

clean-all: clean-samples clean

###########################################################################
## Default targets (called for plain `make`)

default: $(DEFAULT_TARGETS)
