cmake_minimum_required(VERSION 3.20)
project(Summer25_Research LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT FLIPDIST_IPO_SUPPORTED OUTPUT FLIPDIST_IPO_ERROR)

option(FLIPDIST_REUSE_PROFILE_COUNTERS
  "Enable profile-only unique/repeated reuse counters in recursive hot paths"
  OFF)

set(FLIPDIST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/flipdist)

# bf_bst: driver + helpers all in bf_bst.cpp
add_executable(bf_bst
  ${FLIPDIST_SRC_DIR}/bf_bst.cpp
)
target_include_directories(bf_bst PRIVATE ${FLIPDIST_SRC_DIR})

# flipdist: main.cpp supplies main(), bf_bst.cpp supplies helpers
add_executable(flipdist
  ${FLIPDIST_SRC_DIR}/main.cpp
  ${FLIPDIST_SRC_DIR}/bf_bst.cpp
  ${FLIPDIST_SRC_DIR}/memoization.cpp
  ${FLIPDIST_SRC_DIR}/helpers.cpp
  ${FLIPDIST_SRC_DIR}/algorithm.cpp
  ${FLIPDIST_SRC_DIR}/testing_cli.cpp
)
# Skip the bf_bst.cpp main() when building flipdist
target_compile_definitions(flipdist PRIVATE FLIPDIST_BINARY)
if(FLIPDIST_REUSE_PROFILE_COUNTERS)
  target_compile_definitions(flipdist PRIVATE FLIPDIST_REUSE_PROFILE_COUNTERS=1)
endif()
target_include_directories(flipdist PRIVATE
  ${FLIPDIST_SRC_DIR}
)
if(FLIPDIST_IPO_SUPPORTED)
  set_property(TARGET flipdist PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

# “make compile” now correctly depends on both targets
add_custom_target(compile
  DEPENDS bf_bst flipdist
)

include(CTest)
if(BUILD_TESTING)
  add_test(
    NAME flipdist_smoke
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke.sh --build-dir ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
