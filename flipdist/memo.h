// Memoization keys and caches for FlipDist
#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "../rotation_tree.h"
#include "types.h"

// Serialize tree structure into preorder/inorder vectors and a canonical signature.
std::pair<std::vector<int>, std::vector<int>> serializeTree(const VectorRangeTreeMap &T);
const std::string& treeSignature(const VectorRangeTreeMap &T);
const std::string& treeShapeSignature(const VectorRangeTreeMap &T);

struct FlipMemoKey {
    std::string start;
    std::string target;
    int k;

    bool operator==(const FlipMemoKey &other) const noexcept {
        return k == other.k && start == other.start && target == other.target;
    }
};

struct FlipMemoKeyHash {
    std::size_t operator()(const FlipMemoKey &key) const noexcept;
};

extern std::unordered_map<FlipMemoKey, bool, FlipMemoKeyHash> g_flipDistMemo;
extern std::unordered_map<std::string, bool> g_treeDistSMemo;
extern std::unordered_map<std::string, bool> g_treeDistIMemo;

// Monotone bounds cache for TreeDistS: if a subproblem is solvable within k,
// it is solvable for any larger budget; if it is not solvable within k, it is
// not solvable for any smaller budget.
struct BudgetBounds
{
    int min_true = 0x3fffffff; // smallest k known to succeed
    int max_false = -1;        // largest k known to fail
};

extern std::unordered_map<std::string, BudgetBounds> g_treeDistSBounds;
extern std::unordered_map<std::string, BudgetBounds> g_flipDistBounds;

struct PartitionCacheEntry
{
    VectorRangeTreeMap startInside;
    VectorRangeTreeMap startOutside;
    VectorRangeTreeMap targetInside;
    VectorRangeTreeMap targetOutside;
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S_inside;
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S_outside;
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> conflicts_inside;
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> conflicts_outside;
    int edgesInside = 0;
    int edgesOutside = 0;
    std::unordered_map<int, bool> resultMemo;
};

extern std::unordered_set<std::string> g_partitionFailureMemo;
extern std::unordered_map<std::string, PartitionCacheEntry> g_partitionSuccessMemo;

// Free-edge cache key and store.
std::string makeFreeEdgeCacheKey(const VectorRangeTreeMap &start,
                                 const VectorRangeTreeMap &target);
extern std::unordered_map<std::string, std::vector<std::pair<int,int>>> g_freeEdgeCache;
