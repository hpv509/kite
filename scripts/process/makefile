TARGET    := process
PREFIX  := build
CXX     := g++
OBJDIR  := $(PREFIX)/.obj

EIGEN_INC := -I/usr/local/include/eigen3

HDF5_LIBS   := -lhdf5_cpp -lhdf5

CPPFLAGS := $(EIGEN_INC)

CXXFLAGS := -std=c++20 -Ofast
LDFLAGS  := $(HDF5_LIBS)

SRC := $(wildcard *.cpp)

OBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC))

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(PREFIX)
	$(CXX) $(LDFLAGS) $^ -o $(PREFIX)/$@

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(OBJ:.o=.d)

clean:
	rm -rf $(PREFIX)

.PHONY: all clean
