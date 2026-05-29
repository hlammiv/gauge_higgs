# Build for the arbitrary-D SU(N) gauge + arbitrary-irrep Higgs HMC code.
# Header-only library under src/; this builds the driver and the test binaries.
# MPI is not yet used; CXX defaults to g++. Override with `make CXX=mpicxx` once
# MPI domain decomposition lands. CPU-only (no CUDA/cmake on this machine).

CXX      ?= g++
STD      ?= -std=c++20
OPT      ?= -O3 -march=native -funroll-loops
WARN     ?= -Wall -Wextra -Wno-unused-parameter
OMP      ?= -fopenmp
INC      := -Isrc
CXXFLAGS ?= $(STD) $(OPT) $(WARN) $(OMP) $(INC)

# Default lattice config for the driver (override: make NDIM=3 NCOL=3).
NDIM ?= 4
NCOL ?= 2
DEFS := -DNDIM=$(NDIM) -DNCOL=$(NCOL)

BUILD      := build
HEADERS    := $(shell find src -name '*.hpp')
# screening.cpp is a D=3-only demo driver (static_assert NDIM==3); build on demand with
# `make NDIM=3 NCOL=2 build/screening`, not in the default (NDIM=4) `all`.
DRIVERSRC  := $(filter-out src/screening.cpp,$(wildcard src/*.cpp))
DRIVERBINS := $(patsubst src/%.cpp,$(BUILD)/%,$(DRIVERSRC))
TESTSRC    := $(wildcard test/*.cpp)
TESTBINS   := $(patsubst test/%.cpp,$(BUILD)/%,$(TESTSRC))

.PHONY: all test clean
all: $(DRIVERBINS) $(TESTBINS)

# Config stamp: its name encodes NDIM/NCOL, so changing them forces a driver rebuild
# (compile-time config can't be picked up by timestamps alone).
CONFIG_STAMP := $(BUILD)/.config-$(NDIM)-$(NCOL)
$(CONFIG_STAMP): | $(BUILD)
	@rm -f $(BUILD)/.config-* && touch $@

# Drivers use the compile-time lattice config (DEFS); tests instantiate templates directly.
$(BUILD)/%: src/%.cpp $(HEADERS) $(CONFIG_STAMP) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(DEFS) -o $@ $<

# MPI domain-decomposed drivers (separate target; needs mpicxx + -DUSE_MPI).
.PHONY: mpi
mpi: $(BUILD)/gh_hmc_mpi $(BUILD)/gh_higgs_mpi $(BUILD)/gh_higgs_multi_mpi
$(BUILD)/gh_hmc_mpi: src/mpi/gh_hmc_mpi.cpp $(HEADERS) $(CONFIG_STAMP) | $(BUILD)
	mpicxx $(STD) $(OPT) $(WARN) $(OMP) $(INC) $(DEFS) -DUSE_MPI -o $@ $<
$(BUILD)/gh_higgs_mpi: src/mpi/gh_higgs_mpi.cpp $(HEADERS) $(CONFIG_STAMP) | $(BUILD)
	mpicxx $(STD) $(OPT) $(WARN) $(OMP) $(INC) $(DEFS) -DUSE_MPI -o $@ $<
$(BUILD)/gh_higgs_multi_mpi: src/mpi/gh_higgs_multi_mpi.cpp $(HEADERS) $(CONFIG_STAMP) | $(BUILD)
	mpicxx $(STD) $(OPT) $(WARN) $(OMP) $(INC) $(DEFS) -DUSE_MPI -o $@ $<

$(BUILD)/%: test/%.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUILD):
	mkdir -p $(BUILD)

# Build and run every test binary; stop on first failure.
test: $(TESTBINS)
	@rc=0; for t in $(TESTBINS); do \
	  echo "==== $$t ===="; $$t || rc=1; done; \
	exit $$rc

clean:
	rm -rf $(BUILD)
