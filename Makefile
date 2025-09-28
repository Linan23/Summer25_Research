# Makefile (ASan + static lib + optional executables)

CXX      := clang++
SANFLAGS := -fsanitize=address -fno-omit-frame-pointer
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O1 -g $(SANFLAGS)
LDFLAGS  := $(SANFLAGS)

BUILD    := build/asan

# All cpp files
ALL_SRCS  := $(wildcard *.cpp)

# Mains are optional; only link if present
MAIN_SRCS := $(filter FlipDist.cpp test.cpp,$(ALL_SRCS))

# Library sources = everything except mains
LIB_SRCS  := $(filter-out $(MAIN_SRCS),$(ALL_SRCS))
LIB_OBJS  := $(patsubst %.cpp,$(BUILD)/%.o,$(LIB_SRCS))

# Optional main objects
FLIPDIST_MAIN := $(BUILD)/FlipDist.o
TEST_MAIN     := $(BUILD)/test.o

# Default: always build the lib; also build executables if mains exist
.PHONY: all
all: lib \
     $(if $(wildcard FlipDist.cpp),flipdist_asan,) \
     $(if $(wildcard test.cpp),test_asan,)

# Static library with all non-main objects
.PHONY: lib
lib: $(BUILD)/libtreedist.a

$(BUILD)/libtreedist.a: $(LIB_OBJS) | $(BUILD)
	@mkdir -p $(BUILD)
	ar rcs $@ $^

# Executables (only built if their sources exist)
flipdist_asan: $(FLIPDIST_MAIN) $(BUILD)/libtreedist.a
	$(CXX) $^ -o $@ $(LDFLAGS)

test_asan: $(TEST_MAIN) $(BUILD)/libtreedist.a
	$(CXX) $^ -o $@ $(LDFLAGS)

# Generic compile rule
$(BUILD)/%.o: %.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Ensure build dir exists
$(BUILD):
	mkdir -p $(BUILD)

# Clean targets
.PHONY: clean clobber
clean:
	rm -rf $(BUILD)/*.o $(BUILD)/*.d $(BUILD)/libtreedist.a

clobber: clean
	rm -rf build flipdist_asan test_asan

# Auto-include deps
-include $(LIB_OBJS:.o=.d) \
         $(if $(wildcard FlipDist.cpp),$(FLIPDIST_MAIN:.o=.d)) \
         $(if $(wildcard test.cpp),$(TEST_MAIN:.o=.d))
