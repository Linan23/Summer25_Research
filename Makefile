# Makefile (ASan + static lib + optional executables)

CXX      := clang++
SANFLAGS := -fsanitize=address -fno-omit-frame-pointer
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O1 -g $(SANFLAGS)
LDFLAGS  := $(SANFLAGS)

BUILD    := build/asan
REL_BUILD := build/release
REL_CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O3 -DNDEBUG
REL_LDFLAGS :=

# All cpp files
ALL_SRCS  := $(wildcard *.cpp) $(wildcard flipdist/*.cpp)

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

# Release build (no sanitizers, optimized).
.PHONY: release
release: release_lib \
	 $(if $(wildcard FlipDist.cpp),flipdist_fast,) \
	 $(if $(wildcard test.cpp),test_fast,)

# Static library with all non-main objects
.PHONY: lib
lib: $(BUILD)/libtreedist.a

$(BUILD)/libtreedist.a: $(LIB_OBJS) | $(BUILD)
	@mkdir -p $(BUILD)
	ar rcs $@ $^

# Release static library
.PHONY: release_lib
release_lib: $(REL_BUILD)/libtreedist.a

REL_LIB_OBJS := $(patsubst %.cpp,$(REL_BUILD)/%.o,$(LIB_SRCS))
$(REL_BUILD)/libtreedist.a: $(REL_LIB_OBJS) | $(REL_BUILD)
	@mkdir -p $(REL_BUILD)
	ar rcs $@ $^

# Executables (only built if their sources exist)
flipdist_asan: $(FLIPDIST_MAIN) $(BUILD)/libtreedist.a
	$(CXX) $^ -o $@ $(LDFLAGS)

test_asan: $(TEST_MAIN) $(BUILD)/libtreedist.a
	$(CXX) $^ -o $@ $(LDFLAGS)

REL_FLIPDIST_MAIN := $(REL_BUILD)/FlipDist.o
REL_TEST_MAIN     := $(REL_BUILD)/test.o

flipdist_fast: $(REL_FLIPDIST_MAIN) $(REL_BUILD)/libtreedist.a
	$(CXX) $^ -o $@ $(REL_LDFLAGS)

test_fast: $(REL_TEST_MAIN) $(REL_BUILD)/libtreedist.a
	$(CXX) $^ -o $@ $(REL_LDFLAGS)

# Generic compile rule
$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Release compile rule
$(REL_BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(REL_CXXFLAGS) -MMD -MP -c $< -o $@

# Ensure build dir exists
$(BUILD):
	mkdir -p $(BUILD)

$(REL_BUILD):
	mkdir -p $(REL_BUILD)

# Clean targets
.PHONY: clean clobber
clean:
	rm -rf $(BUILD)/*.o $(BUILD)/*.d $(BUILD)/libtreedist.a

clobber: clean
	rm -rf build flipdist_asan test_asan flipdist_fast test_fast

# Auto-include deps
-include $(LIB_OBJS:.o=.d) \
         $(if $(wildcard FlipDist.cpp),$(FLIPDIST_MAIN:.o=.d)) \
         $(if $(wildcard test.cpp),$(TEST_MAIN:.o=.d))

-include $(REL_LIB_OBJS:.o=.d) \
	 $(if $(wildcard FlipDist.cpp),$(REL_FLIPDIST_MAIN:.o=.d)) \
	 $(if $(wildcard test.cpp),$(REL_TEST_MAIN:.o=.d))
