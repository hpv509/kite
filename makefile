PREFIX  := build
CXX     := g++
LINK    := h5c++
OBJDIR  := $(PREFIX)/.obj

# Eigen-Dir
EIGEN_INC := -I/usr/local/include/eigen3

# HDF5-Dir
HDF5_FLAGS := $(shell h5c++ -show | tr ' ' '\n' | grep '^-I' | xargs)

CPPFLAGS := $(EIGEN_INC) \
	-Itools/Src \
	-ISrc \
	-ISrc/Tools \
	-ISrc/Lattice \
	-ISrc/Hamiltonian \
	-ISrc/Vector \
	-ISrc/Simulation \
	-DCOMPILE_WAVEPACKET=0 \
	$(HDF5_FLAGS)

CXXFLAGS := -std=c++17 -O3 -fopenmp
LDFLAGS  := -fopenmp

SRC_KITEX := $(shell find Src -name '*.cpp')
SRC_TOOLS := $(shell find tools/Src -name '*.cpp')

OBJ_KITEX := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_KITEX))
OBJ_TOOLS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_TOOLS))

all: KITEx KITE-tools

KITEx: $(OBJ_KITEX)
	@mkdir -p $(PREFIX)
	$(LINK) $(LDFLAGS) $^ -o $(PREFIX)/$@

KITE-tools: $(OBJ_TOOLS)
	@mkdir -p $(PREFIX)
	$(LINK) $(LDFLAGS) $^ -o $(PREFIX)/$@

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(OBJ_KITEX:.o=.d) $(OBJ_TOOLS:.o=.d)

clean:
	rm -rf $(PREFIX)

.PHONY: all clean
