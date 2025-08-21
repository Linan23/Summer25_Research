cmake_minimum_required(VERSION 3.30)
project(Summer25_Research LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# A-tree: driver + helpers all in A_tree.cpp
add_executable(A_tree
  A_tree.cpp
)

# B-tree: FlipDist.cpp supplies main(), A_tree.cpp supplies helpers
add_executable(B_tree
  FlipDist.cpp
  A_tree.cpp
)
# Skip the A_tree.cpp main() when building B_tree
target_compile_definitions(B_tree PRIVATE FLIPDIST_BINARY)

# “make compile” now correctly depends on both targets
add_custom_target(compile
  DEPENDS A_tree B_tree
)
