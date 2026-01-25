// Branching helpers for FlipDist (pair filtering, partitioning)
#pragma once

#include <vector>
#include <utility>
#include "../rotation_tree.h"
#include "types.h"

// Validate and append a diagonal pair.
bool appendValidPair(std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &dest,
                     const DiagonalEdge &a,
                     const DiagonalEdge &b);

// Remove duplicate diagonal pairs (order-insensitive within the pair and list).
void deduplicateDiagonalPairs(std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &pairs,
                              const char *tag);

// Lightweight profiling hook.
void recordBranchSample(const char *tag, size_t before, size_t after);

// Partition S between two disjoint subtrees.
std::pair<std::vector<std::pair<DiagonalEdge, DiagonalEdge>>,
          std::vector<std::pair<DiagonalEdge, DiagonalEdge>>>
partitionS(const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S,
           const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2);

