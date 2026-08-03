#pragma once

#include "bf_bst.h"

#include <utility>
#include <vector>

// Purpose: Tests whether two trees are within the supplied rotation budget.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); k (int).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool FlipDistTree(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k);

// Purpose: Tests a branch defined by an independent set of rotations.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); k (int); I (const std::vector<std::pair<int, int>> &).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool TreeDistI(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k,
               const std::vector<std::pair<int, int>> &I);

// Purpose: Tests the remaining budget using free-edge decomposition or S-pair branching.
// Inputs: T_init (const VectorRangeTreeMap &); T_end (const VectorRangeTreeMap &); k (int); S (const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool TreeDistS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
               const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);

// Purpose: Finds the minimum feasible rotation budget up to the supplied cap.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); max_k (int).
// Output: The minimum feasible budget, or -1 when no solution is found through max_k.
int FlipDistMinK(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int max_k);

// Purpose: Resets algorithm caches.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void resetAlgorithmCaches();

// Purpose: Resets algorithm profile sets.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void resetAlgorithmProfileSets();
