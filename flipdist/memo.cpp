// Memoization helpers for FlipDist
#include "memo.h"
#include "utils.h"
#include "profile.h"
#include <algorithm>

std::unordered_map<FlipMemoKey, bool, FlipMemoKeyHash> g_flipDistMemo;
std::unordered_map<std::string, bool> g_treeDistSMemo;
std::unordered_map<std::string, bool> g_treeDistIMemo;
std::unordered_map<std::string, BudgetBounds> g_treeDistSBounds;
std::unordered_map<std::string, BudgetBounds> g_flipDistBounds;
std::unordered_set<std::string> g_partitionFailureMemo;
std::unordered_map<std::string, PartitionCacheEntry> g_partitionSuccessMemo;
std::unordered_map<std::string, std::vector<std::pair<int,int>>> g_freeEdgeCache;

std::size_t FlipMemoKeyHash::operator()(const FlipMemoKey &key) const noexcept {
    std::size_t h1 = std::hash<std::string>()(key.start);
    std::size_t h2 = std::hash<std::string>()(key.target);
    std::size_t h3 = std::hash<int>()(key.k);
    std::size_t seed = h1;
    seed ^= h2 + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
    return seed;
}

std::pair<std::vector<int>, std::vector<int>> serializeTree(const VectorRangeTreeMap &T)
{
    std::vector<int> preorder;
    std::vector<int> inorder;

    std::function<void(int)> dfsPre = [&](int node)
    {
        if (node < 0 || !T.isOriginal(node)) return;
        preorder.push_back(node);
        dfsPre(T.getLeftChild(node));
        dfsPre(T.getRightChild(node));
    };

    std::function<void(int)> dfsIn = [&](int node)
    {
        if (node < 0 || !T.isOriginal(node)) return;
        dfsIn(T.getLeftChild(node));
        inorder.push_back(node);
        dfsIn(T.getRightChild(node));
    };

    if (T.root >= 0 && T.isOriginal(T.root))
    {
        dfsPre(T.root);
        dfsIn(T.root);
    }

    return {preorder, inorder};
}

const std::string& treeSignature(const VectorRangeTreeMap &T)
{
    // Prefer the tree's cached signature (invalidated on rotation/build) to avoid
    // repeatedly serialising the same structure during memoised recursion.
    flipdist::profile::ScopedTimer timer(&flipdist::profile::g.ns_treeSignature);
    if (flipdist::profile::enabled())
        ++flipdist::profile::g.calls_treeSignature;
    return T.signature();
}

const std::string& treeShapeSignature(const VectorRangeTreeMap &T)
{
    // Shape signature is canonical under relabeling; use it for memo keys that
    // do not store concrete node IDs (e.g., boolean subproblem results).
    return T.shapeSignature();
}

std::string makeFreeEdgeCacheKey(const VectorRangeTreeMap &start,
                                 const VectorRangeTreeMap &target)
{
    const auto &a = treeSignature(start);
    const auto &b = treeSignature(target);
    std::string key;
    key.reserve(a.size() + b.size() + 2);
    key.append(a);
    key.append("||");
    key.append(b);
    return key;
}
