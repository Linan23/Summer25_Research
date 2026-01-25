// Conflict detection and edge utilities for FlipDist
#include "conflicts.h"
#include "treedist.h"
#include "profile.h"
#include "memo.h"
#include "../rotation_tree.h"
#include "partners.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace {
bool findEdgeOrientationLocal(const VectorRangeTreeMap &T,
                              int a,
                              int b,
                              int &parent_out,
                              int &child_out)
{
    if (!T.isOriginal(a) || !T.isOriginal(b)) return false;
    if (T.getLeftChild(a) == b || T.getRightChild(a) == b) {
        parent_out = a;
        child_out = b;
        return true;
    }
    if (T.getLeftChild(b) == a || T.getRightChild(b) == a) {
        parent_out = b;
        child_out = a;
        return true;
    }
    return false;
}
} // namespace

static inline std::pair<int,int> makeUndirectedPair(int a, int b)
{
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

std::vector<std::pair<int, int>> getInternalEdges(const VectorRangeTreeMap &T)
{
    std::vector<std::pair<int, int>> edges;
    try
    {
        if (T.original_nodes.empty() || T.root < 0)
        {
            return edges;
        }

        std::function<void(int)> dfs = [&](int node)
        {
            if (node < 0 || !T.isOriginal(node))
                return;
            try
            {
                int left = T.getLeftChild(node);
                int right = T.getRightChild(node);
                if (left >= 0 && T.isOriginal(left))
                {
                    edges.emplace_back(node, left);
                    dfs(left);
                }
                if (right >= 0 && T.isOriginal(right))
                {
                    edges.emplace_back(node, right);
                    dfs(right);
                }
            }
            catch (...)
            {
                return;
            }
        };

        if (T.isOriginal(T.root))
        {
            dfs(T.root);
        }
    }
    catch (...)
    {
        edges.clear();
    }
    return edges;
}

std::vector<std::pair<int, int>> getIncidentEdges(const VectorRangeTreeMap &T, int node)
{
    std::vector<std::pair<int, int>> incident;
    try
    {
        if (!T.isOriginal(node))
            return incident;

        int left = T.getLeftChild(node);
        int right = T.getRightChild(node);
        if (left >= 0 && T.isOriginal(left))
        {
            incident.emplace_back(node, left);
        }
        if (right >= 0 && T.isOriginal(right))
        {
            incident.emplace_back(node, right);
        }

        int parent = T.getParent(node);
        if (parent >= 0 && T.isOriginal(parent))
        {
            incident.emplace_back(parent, node);
        }
    }
    catch (...)
    {
        incident.clear();
    }
    return incident;
}

int lowerBoundEdgeDifference(const VectorRangeTreeMap &A,
                             const VectorRangeTreeMap &B)
{
    // Lower bound based on subtree ranges (intervals in inorder space) rather
    // than parent→child labels. To keep this hot-path cheap, represent each
    // tree's range-set as a compact bitmask over all (L,R) intervals in
    // {0..m} where m = #nodes (endpoints = m+1)
    const int m = static_cast<int>(A.original_inorder.size());
    if (m <= 0)
        return 0;
    if (static_cast<int>(B.original_inorder.size()) != m)
        return 0;

// Cache lower-bound computations keyed by shape signatures (symmetric key)
static std::unordered_map<std::string, int> lbCache;
    static std::size_t lbCap = 0;
    if (lbCap == 0)
    {
        if (const char *env = std::getenv("FLIPDIST_LB_CACHE_CAP"))
        {
            long long v = std::atoll(env);
            if (v > 0)
                lbCap = static_cast<std::size_t>(v);
        }
        if (lbCap == 0)
            lbCap = 200000;
    }

    const auto &sigA = treeShapeSignature(A);
    const auto &sigB = treeShapeSignature(B);
    std::string key;
    key.reserve(sigA.size() + sigB.size() + 2);
    if (sigA <= sigB)
    {
        key.append(sigA);
        key.append("||");
        key.append(sigB);
    }
    else
    {
        key.append(sigB);
        key.append("||");
        key.append(sigA);
    }

    if (auto it = lbCache.find(key); it != lbCache.end())
        return it->second;

    const int intervalCount = (m + 1) * m / 2;
    const int wordCount = (intervalCount + 63) / 64;

    auto intervalIndex = [&](int L, int R) -> int {
        // 0 <= L < R <= m.
        const int prefix = L * m - (L * (L - 1)) / 2;
        return prefix + (R - L - 1);
    };

    int diff = 0;
    if (wordCount <= 6)
    {
        auto buildMask = [&](const VectorRangeTreeMap &T,
                             std::array<std::uint64_t, 6> &mask) {
            mask.fill(0);
            for (int node : T.original_inorder)
            {
                auto r = T.getRange(node);
                if (r.first < 0 || r.second <= r.first)
                    continue;
                if (r.first > m || r.second > m)
                    continue;
                int idx = intervalIndex(r.first, r.second);
                int word = idx >> 6;
                int bit = idx & 63;
                mask[static_cast<size_t>(word)] |= (1ull << bit);
            }
        };

        std::array<std::uint64_t, 6> maskA{};
        std::array<std::uint64_t, 6> maskB{};
        buildMask(A, maskA);
        buildMask(B, maskB);
        for (int i = 0; i < wordCount; ++i)
        {
            diff += __builtin_popcountll(maskA[static_cast<size_t>(i)] ^ maskB[static_cast<size_t>(i)]);
        }
    }
    else
    {
        auto buildMask = [&](const VectorRangeTreeMap &T,
                             std::vector<std::uint64_t> &mask) {
            std::fill(mask.begin(), mask.end(), 0);
            for (int node : T.original_inorder)
            {
                auto r = T.getRange(node);
                if (r.first < 0 || r.second <= r.first)
                    continue;
                if (r.first > m || r.second > m)
                    continue;
                int idx = intervalIndex(r.first, r.second);
                int word = idx >> 6;
                int bit = idx & 63;
                mask[static_cast<size_t>(word)] |= (1ull << bit);
            }
        };

        std::vector<std::uint64_t> maskA(static_cast<size_t>(wordCount), 0);
        std::vector<std::uint64_t> maskB(static_cast<size_t>(wordCount), 0);
        buildMask(A, maskA);
        buildMask(B, maskB);
        for (int i = 0; i < wordCount; ++i)
        {
            diff += __builtin_popcountll(maskA[static_cast<size_t>(i)] ^ maskB[static_cast<size_t>(i)]);
        }
    }

    // In this subtree-range (interval) model, a single rotation replaces exactly
    // one interval in the range-set (the old child's interval) with a new one
    // (the demoted parent's interval after rotation), so the symmetric
    // difference can shrink by at most 2 per rotation.
    int result = diff / 2;
    if (lbCache.size() > lbCap)
        lbCache.clear();
    lbCache.emplace(std::move(key), result);
    return result;
}

int countInternalEdges(const VectorRangeTreeMap &T)
{
    return static_cast<int>(getInternalEdges(T).size());
}

std::vector<DiagonalEdge> collectConflictingEdges(
        const VectorRangeTreeMap &start,
        const VectorRangeTreeMap &target)
{
    flipdist::profile::ScopedTimer timer(&flipdist::profile::g.ns_collectConflicts);
    if (flipdist::profile::enabled())
        ++flipdist::profile::g.calls_collectConflicts;

    auto startIndex  = buildEndpointIndex(start);
    auto targetIndex = buildEndpointIndex(target);

    std::vector<DiagonalEdge> conflicts;
    conflicts.reserve(startIndex.size());

    auto makeOriented = [&](int parent, int child, const std::pair<int,int>& diag) -> DiagonalEdge {
        return DiagonalEdge{diag, parent, child};
    };

    for (const auto &entry : startIndex)
    {
        const auto &diag = entry.first;
        if (targetIndex.count(diag))
            continue;

        int node = entry.second;
        if (!start.isOriginal(node))
            continue;

        int parent = start.getParent(node);
        if (parent != VectorRangeTreeMap::NO_PARENT && start.isOriginal(parent))
        {
            conflicts.push_back(makeOriented(parent, node, diag));
            continue;
        }

        int left = start.getLeftChild(node);
        if (left != VectorRangeTreeMap::NO_CHILD && start.isOriginal(left))
        {
            conflicts.push_back(makeOriented(node, left, diag));
            continue;
        }

        int right = start.getRightChild(node);
        if (right != VectorRangeTreeMap::NO_CHILD && start.isOriginal(right))
        {
            conflicts.push_back(makeOriented(node, right, diag));
            continue;
        }
    }

    if (conflicts.empty())
    {
        // Fallback to edge-based difference when diagonals coincide but orientation differs.
        std::unordered_set<std::pair<int,int>, PairHash, PairEq> targetEdgesDirected;
        target.collectEdges(target.root, targetEdgesDirected);
        std::unordered_set<std::pair<int,int>, PairHash, PairEq> targetUndirected;
        std::unordered_map<std::pair<int,int>, std::pair<int,int>, PairHash, PairEq> targetOrient;
        for (const auto &edge : targetEdgesDirected)
        {
            targetUndirected.insert(makeUndirectedPair(edge.first, edge.second));
            int p, c;
            if (findEdgeOrientationLocal(target, edge.first, edge.second, p, c))
                targetOrient[makeUndirectedPair(edge.first, edge.second)] = std::make_pair(p, c);
        }

        std::unordered_set<std::pair<int,int>, PairHash, PairEq> startEdgesDirected;
        start.collectEdges(start.root, startEdgesDirected);
        for (const auto &edge : startEdgesDirected)
        {
            auto und = makeUndirectedPair(edge.first, edge.second);
            if (!targetUndirected.count(und))
            {
                auto diag = start.diagonalEndpoints(edge.second);
                conflicts.push_back(makeOriented(edge.first, edge.second, diag));
            }
            else
            {
                // Same undirected edge exists, but if orientation differs, treat as conflict.
                int p, c;
                if (findEdgeOrientationLocal(start, edge.first, edge.second, p, c))
                {
                    auto it = targetOrient.find(und);
                    if (it != targetOrient.end())
                    {
                        if (it->second.first != p || it->second.second != c)
                        {
                            auto diag = start.diagonalEndpoints(edge.second);
                            conflicts.push_back(makeOriented(p, c, diag));
                        }
                    }
                }
            }
        }
    }

    if (flipdist::profile::enabled())
        flipdist::profile::g.max_conflicts_seen = std::max<std::uint64_t>(flipdist::profile::g.max_conflicts_seen,
                                                                          static_cast<std::uint64_t>(conflicts.size()));
    return conflicts;
}

std::vector<std::pair<int,int>> buildMaxIndependentSet(
        const VectorRangeTreeMap &start,
        const std::vector<DiagonalEdge> &conflicts)
{
    std::vector<std::pair<int,int>> independent;
    std::unordered_set<int> usedNodes;

    std::vector<DiagonalEdge> sorted = conflicts;
    std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        int widthA = a.diag.second - a.diag.first;
        int widthB = b.diag.second - b.diag.first;
        if (widthA != widthB) return widthA > widthB;
        if (a.parent != b.parent) return a.parent < b.parent;
        return a.child < b.child;
    });

    for (const auto &edge : sorted) {
        int parent = edge.parent;
        int child  = edge.child;
        if (!start.isOriginal(parent) || !start.isOriginal(child))
            continue;

        if (usedNodes.count(parent) || usedNodes.count(child))
            continue;

        independent.emplace_back(parent, child);
        usedNodes.insert(parent);
        usedNodes.insert(child);
    }

    return independent;
}
