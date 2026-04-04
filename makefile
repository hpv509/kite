PREFIX  := build
CXX     := g++
OBJDIR  := $(PREFIX)/.obj

# Eigen-Dir
EIGEN_INC := -I /usr/local/include/eigen3

# HDF5 flags (use the ones from pkg-config by default. If absent use fallback option)
HDF5_CXXFLAGS := $(shell pkg-config --cflags hdf5_cpp hdf5 2>/dev/null)
HDF5_LIBS     := $(shell pkg-config --libs   hdf5_cpp hdf5 2>/dev/null)
ifeq ($(strip $(HDF5_LIBS)),)
  HDF5_LIBS := -lhdf5_cpp -lhdf5
endif

CPPFLAGS := $(EIGEN_INC) \
	-I tools/Src \
	-I Src \
	-I Src/Tools \
	-I Src/Lattice \
	-I Src/Hamiltonian \
	-I Src/Vector \
	-I Src/Simulation \
	-DCOMPILE_WAVEPACKET=0 \
	$(HDF5_CXXFLAGS)

CXXFLAGS := -std=c++17 -O3 -fopenmp
LDFLAGS  := -fopenmp
LDLIBS   := $(HDF5_LIBS)

SRC_KITEX := $(shell find Src -name '*.cpp')
SRC_TOOLS := $(shell find tools/Src -name '*.cpp')

OBJ_KITEX := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_KITEX))
OBJ_TOOLS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_TOOLS))

all: KITEx KITE-tools

KITEx: $(OBJ_KITEX)
	@mkdir -p $(PREFIX)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $(PREFIX)/$@

KITE-tools: $(OBJ_TOOLS)
	@mkdir -p $(PREFIX)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $(PREFIX)/$@

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(OBJ_KITEX:.o=.d) $(OBJ_TOOLS:.o=.d)

clean:
	rm -rf $(PREFIX)

.PHONY: all clean
