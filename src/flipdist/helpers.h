#pragma once

#include "bf_bst.h"

#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

// Purpose: Retrieves internal edges.
// Inputs: T (const VectorRangeTreeMap &).
// Output: The resulting std::vector<std::pair<int, int>> value.
std::vector<std::pair<int, int>> getInternalEdges(const VectorRangeTreeMap &T);

// Purpose: Counts internal edges.
// Inputs: T (const VectorRangeTreeMap &).
// Output: The computed value or status code.
int countInternalEdges(const VectorRangeTreeMap &T);

// Purpose: Counts conflict edges.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &).
// Output: The computed value or status code.
int countConflictEdges(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);

// Purpose: Checks whether two edges share an endpoint.
// Inputs: e1 (const std::pair<int, int> &); e2 (const std::pair<int, int> &).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool areAdjacent(const std::pair<int, int> &e1, const std::pair<int, int> &e2);

// Purpose: Copies a tree while rebuilding metadata needed by the search.
// Inputs: T (const VectorRangeTreeMap &).
// Output: The resulting VectorRangeTreeMap value.
VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap &T);

// Purpose: Reports whether parent child edge.
// Inputs: T (const VectorRangeTreeMap &); parent (int); child (int).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool hasParentChildEdge(const VectorRangeTreeMap &T, int parent, int child);

// Purpose: Resolves an undirected edge to its parent-to-child orientation in a tree.
// Inputs: T (const VectorRangeTreeMap &); edge (const std::pair<int, int> &).
// Output: The resulting std::pair<int, int> value.
std::pair<int, int> orientEdge(const VectorRangeTreeMap &T, const std::pair<int, int> &edge);

// Purpose: Finds free edge.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &).
// Output: The resulting std::pair<bool, std::pair<int, int>> value.
std::pair<bool, std::pair<int, int>> findFreeEdge(const VectorRangeTreeMap &T_init,
                                                  const VectorRangeTreeMap &T_final);

// Purpose: Attempts common edge decomposition.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); k (int); handled (bool &).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool tryCommonEdgeDecomposition(const VectorRangeTreeMap &T_init,
                                const VectorRangeTreeMap &T_final,
                                int k,
                                bool &handled);

// Purpose: Generates all independent subsets.
// Inputs: edges (const std::vector<std::pair<int, int>> &); index (int); current (std::vector<std::pair<int, int>> &); result (std::vector<std::vector<std::pair<int, int>>> &).
// Output: No return value; updates referenced outputs or shared state.
void generateAllIndependentSubsets(const std::vector<std::pair<int, int>> &edges, int index,
                                   std::vector<std::pair<int, int>> &current,
                                   std::vector<std::vector<std::pair<int, int>>> &result);

// Purpose: Retrieves incident edges.
// Inputs: T (const VectorRangeTreeMap &); node (int).
// Output: The resulting std::vector<std::pair<int, int>> value.
std::vector<std::pair<int, int>> getIncidentEdges(const VectorRangeTreeMap &T, int node);

// Purpose: Assigns each S-pair wholly contained in one of two tree partitions.
// Inputs: S: edge pairs to partition; T1 and T2: the two partition trees.
// Output: The S-pairs belonging to T1 and T2, respectively; cross-partition pairs are omitted.
std::pair<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>,
          std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>
partitionS(const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S,
           const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2);

// Purpose: Checks whether two trees produce the same node partition for a child range.
// Inputs: left (const VectorRangeTreeMap&); right (const VectorRangeTreeMap&); child_range (const std::pair<int, int>&).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool partitionNodeSetsMatchByChildRange(const VectorRangeTreeMap& left,
                                        const VectorRangeTreeMap& right,
                                        const std::pair<int, int>& child_range);

// Purpose: Splits S-pairs according to whether their endpoints fall inside or outside a child range.
// Inputs: S: edge pairs to split; T: source tree; child_range: inorder interval defining one partition.
// Output: The S-pairs assigned to the inside and outside partitions, respectively.
std::pair<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>,
          std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>
partitionSByChildRange(const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& S,
                       const VectorRangeTreeMap& T,
                       const std::pair<int, int>& child_range);

// Purpose: Generates independent subsets from s.
// Inputs: S (const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &).
// Output: The resulting std::vector<std::vector<std::pair<int, int>>> value.
std::vector<std::vector<std::pair<int, int>>> generateIndependentSubsetsFromS(
    const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);

struct PrefixPairgenRequestStats {
    int vcount = 0;
    int bucket_size_l = 0;
    int bucket_size_r = 0;
    int emitted_pairs = 0;
};

// Purpose: Appends partner pairs from diagonals.
// Inputs: T (const VectorRangeTreeMap &); diagonals (const std::vector<std::pair<int, int>> &); out_pairs (std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &).
// Output: No return value; updates referenced outputs or shared state.
void appendPartnerPairsFromDiagonals(
    const VectorRangeTreeMap &T,
    const std::vector<std::pair<int, int>> &diagonals,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &out_pairs);

// Purpose: Appends partner pairs from single diagonal profiled.
// Inputs: T (const VectorRangeTreeMap &); diagonal (const std::pair<int, int> &); out_pairs (std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &); stats (PrefixPairgenRequestStats *).
// Output: No return value; updates referenced outputs or shared state.
void appendPartnerPairsFromSingleDiagonalProfiled(
    const VectorRangeTreeMap &T,
    const std::pair<int, int> &diagonal,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &out_pairs,
    PrefixPairgenRequestStats *stats);
