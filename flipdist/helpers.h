#pragma once

#include "bf_bst.h"

#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

std::vector<std::pair<int, int>> getInternalEdges(const VectorRangeTreeMap &T);
int countInternalEdges(const VectorRangeTreeMap &T);
int countConflictEdges(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);
bool areAdjacent(const std::pair<int, int> &e1, const std::pair<int, int> &e2);

VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap &T);
bool hasParentChildEdge(const VectorRangeTreeMap &T, int parent, int child);
std::pair<int, int> orientEdge(const VectorRangeTreeMap &T, const std::pair<int, int> &edge);
std::pair<bool, std::pair<int, int>> findFreeEdge(const VectorRangeTreeMap &T_init,
                                                  const VectorRangeTreeMap &T_final);

bool tryCommonEdgeDecomposition(const VectorRangeTreeMap &T_init,
                                const VectorRangeTreeMap &T_final,
                                int k,
                                bool &handled);

void generateAllIndependentSubsets(const std::vector<std::pair<int, int>> &edges, int index,
                                   std::vector<std::pair<int, int>> &current,
                                   std::vector<std::vector<std::pair<int, int>>> &result);

std::vector<std::pair<int, int>> getIncidentEdges(const VectorRangeTreeMap &T, int node);

std::pair<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>,
          std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>
partitionS(const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S,
           const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2);

std::vector<std::vector<std::pair<int, int>>> generateIndependentSubsetsFromS(
    const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);

void appendPartnerPairsFromDiagonals(
    const VectorRangeTreeMap &T,
    const std::vector<std::pair<int, int>> &diagonals,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &out_pairs);
