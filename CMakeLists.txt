cmake_minimum_required(VERSION 3.30)
project(Summer25_Research LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build A-tree binary
add_executable(A_tree
  A_tree.cpp
)

# Build B-tree binary
add_executable(B-tree
  B-tree.cpp
)

# “make compile” builds *both*
add_custom_target(compile
  DEPENDS A_tree B-tree
)
