cmake_minimum_required(VERSION 3.30)
project(Summer25_Research LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# bf_bst: driver + helpers all in bf_bst.cpp
add_executable(bf_bst
  bf_bst.cpp
)

# flipdist: FlipDist.cpp supplies main(), bf_bst.cpp supplies helpers
add_executable(flipdist
  FlipDist.cpp
  bf_bst.cpp
)
# Skip the bf_bst.cpp main() when building flipdist
target_compile_definitions(flipdist PRIVATE FLIPDIST_BINARY)

# “make compile” now correctly depends on both targets
add_custom_target(compile
  DEPENDS bf_bst flipdist
)
