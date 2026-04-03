include config.mk
OBJDIR := build/.obj

SRC_KITEX  := $(shell find Src -name '*.cpp')
SRC_TOOLS  := $(shell find tools/Src -name '*.cpp')
OBJ_KITEX := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_KITEX))
OBJ_TOOLS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC_TOOLS))

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

all: KITEx KITE-tools

KITEx: $(OBJ_KITEX)
	@mkdir -p build
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o build/KITEx

KITE-tools: $(OBJ_TOOLS)
	@mkdir -p build
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o build/KITE-tools

build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(OBJ_KITEX:.o=.d) $(OBJ_TOOLS:.o=.d)

clean:
	rm -rf build

.PHONY: all clean
