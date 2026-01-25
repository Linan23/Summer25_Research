// Core flip distance algorithm entry points
#pragma once

#include <unordered_set>
#include <vector>
#include "../rotation_tree.h"
#include "types.h"

namespace flipdist {

int FlipDistMinK(const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2, int k_max);
bool FlipDistTree(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k);
bool TreeDistI(const VectorRangeTreeMap &T_init,
               const VectorRangeTreeMap &T_final,
               int k,
               const std::vector<std::pair<int, int>> &I,
               bool allow_independent_retry = true,
               const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflictSet = nullptr);
bool TreeDistS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
               const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S,
               bool allow_independent_retry = true,
               const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflictSet = nullptr);
bool hasParentChildEdge(const VectorRangeTreeMap &T, int parent, int child);

} // namespace flipdist
