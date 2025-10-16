cmake_minimum_required(VERSION 3.16)
project(Summer25_Research LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Common library with all core code (no main)
add_library(flipcore STATIC
  tree_core.cpp
  tree_rotate.cpp
  serialize.cpp
  free_edge.cpp
  bfs.cpp
  distance.cpp
  test.cpp         # ok to include; it only provides helpers, not main
)

# So headers like "rotation_tree.h" are found by anything linking to flipcore
target_include_directories(flipcore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


# B_tree: your FlipDist program (has its own main in FlipDist.cpp)
add_executable(B_tree
  FlipDist.cpp
)
target_link_libraries(B_tree PRIVATE flipcore)

add_executable(tests test.cpp)
target_link_libraries(tests PRIVATE flipcore)

# Convenience meta-target
add_custom_target(compile
  DEPENDS A_tree B_tree
)

# warnings
if (MSVC)
  target_compile_options(flipcore PRIVATE /W4)
else()
  target_compile_options(flipcore PRIVATE -Wall -Wextra -Wpedantic)
endif()
