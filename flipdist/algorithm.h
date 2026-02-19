#pragma once

#include "bf_bst.h"

#include <utility>
#include <vector>

bool FlipDistTree(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k);
bool TreeDistI(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k,
               const std::vector<std::pair<int, int>> &I);
bool TreeDistS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
               const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);
int FlipDistMinK(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int max_k);
