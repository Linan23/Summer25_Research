// Partner generation and diagonal index helpers for FlipDist
#pragma once

#include <vector>
#include <unordered_map>
#include <utility>
#include "../rotation_tree.h"
#include "types.h"

// Map from diagonal endpoints (L,R) to the node realising it.
std::unordered_map<std::pair<int,int>, int, PairHash, PairEq>
buildEndpointIndex(const VectorRangeTreeMap &tree);

// Thin wrapper around buildEndpointIndex.
std::unordered_map<std::pair<int,int>, int, PairHash, PairEq>
buildDiagonalNodeMap(const VectorRangeTreeMap &tree);

// Generate Li–Xia partner wedges for the given trees.
std::vector<std::pair<DiagonalEdge, DiagonalEdge>>
buildPartnerPairs(const VectorRangeTreeMap &start,
                  const VectorRangeTreeMap &target);
