// Core flip distance algorithm implementation (TreeDistI/TreeDistS/FlipDistTree)
#include "treedist.h"
#include "partners.h"
#include "branch_utils.h"
#include "conflicts.h"
#include "memo.h"
#include "utils.h"
#include "types.h"
#include "profile.h"
#include "../rotation_tree.h"

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <deque>
#include <functional>
#include <optional>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>

extern const bool DEBUG;
void debugPrint(const std::string &msg);
// BFS helper lives in FlipDist.cpp
extern int MinRotationsBFS(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap);

namespace flipdist {
namespace {

static bool earContractionEnabledFor(int lb_conflicts)
{
    static int mode = -2; // -2=uninitialized, -1=auto, 0=off, 1=on
    if (mode == -2)
    {
        if (const char *env = std::getenv("FLIPDIST_EAR_CONTRACTION"))
        {
            int v = std::atoi(env);
            mode = (v > 0) ? 1 : 0;
        }
        else
        {
            mode = -1; // auto
        }
    }

    if (mode == 0)
        return false;
    if (mode == 1)
        return true;

    static int threshold = -1;
    if (threshold == -1)
    {
        if (const char *env = std::getenv("FLIPDIST_EAR_CONTRACTION_LB"))
            threshold = std::max(0, std::atoi(env));
        else
            threshold = 6;
    }

    return lb_conflicts >= threshold;
}
[[maybe_unused]] constexpr bool PROFILE = false;

[[maybe_unused]] static inline bool isInternalDiagonal(const std::pair<int,int> &diag)
{
    return diag.second - diag.first > 1;
}

static inline bool shouldUseInternalBFS(const VectorRangeTreeMap &T, int k)
{
    if (k <= 0)
        return false;
    // Only allow BFS on very small subproblems; otherwise it dominates runtime.
    //
    // NOTE: For >=9 internal edges, even a capped BFS is often slower than the
    // Li–Xia recursion, and it gets invoked many times in hard instances.
    // Keeping this at 8 leverages the precomputed distance table (n<=9 nodes)
    // in `smallRotationDistance`, which is essentially O(1) per query.
    return countInternalEdges(T) <= 8;
}

static int profiledMinRotationsBFS(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap)
{
    struct BfsCacheEntry
    {
        int known_dist = -1;      // exact distance if known (>=0)
        int max_failed_cap = -1;  // largest cap we've proven insufficient
    };

    // Cache bounded BFS outcomes for small subtrees; TreeDistS can revisit the
    // same small pair with varying caps.
    static std::unordered_map<std::string, BfsCacheEntry> bfsMemo;
    profile::ScopedTimer timer(&profile::g.ns_minRotationsBFS_internal);
    if (profile::enabled())
        ++profile::g.calls_minRotationsBFS_internal;

    // Cache small-subtree BFS results; these get called repeatedly inside recursion.
    const int edgesA = countInternalEdges(A);
    const int edgesB = countInternalEdges(B);
    if (edgesA <= 8 && edgesB <= 8 && cap <= 21)
    {
        const auto &sigA = treeShapeSignature(A);
        const auto &sigB = treeShapeSignature(B);

        std::string key;
        key.reserve(sigA.size() + sigB.size() + 2);
        key.append(sigA);
        key.append("||");
        key.append(sigB);

        if (auto it = bfsMemo.find(key); it != bfsMemo.end())
        {
            const auto &entry = it->second;
            if (entry.known_dist >= 0)
            {
                if (profile::enabled())
                    ++profile::g.memo_hits_minRotationsBFS_internal;
                return (entry.known_dist <= cap) ? entry.known_dist : -1;
            }
            if (cap <= entry.max_failed_cap)
            {
                if (profile::enabled())
                    ++profile::g.memo_hits_minRotationsBFS_internal;
                return -1;
            }
        }

        int dist = MinRotationsBFS(A, B, cap);
        if (bfsMemo.size() > 50'000)
            bfsMemo.clear();

        auto &entry = bfsMemo[key];
        if (dist >= 0)
        {
            entry.known_dist = dist;
        }
        else
        {
            entry.known_dist = -1;
            entry.max_failed_cap = std::max(entry.max_failed_cap, cap);
        }

        std::string rev;
        rev.reserve(sigA.size() + sigB.size() + 2);
        rev.append(sigB);
        rev.append("||");
        rev.append(sigA);
        bfsMemo[std::move(rev)] = entry;

        if (profile::enabled())
        {
            if (dist >= 0)
                ++profile::g.successes_minRotationsBFS_internal;
            else
                ++profile::g.failures_minRotationsBFS_internal;
        }
        return dist;
    }

    int dist = MinRotationsBFS(A, B, cap);
    if (profile::enabled())
    {
        if (dist >= 0)
            ++profile::g.successes_minRotationsBFS_internal;
        else
            ++profile::g.failures_minRotationsBFS_internal;
    }
    return dist;
}

// Constructs a DiagonalEdge from an oriented parent->child edge in the given tree.
static DiagonalEdge makeDiagonalEdge(const VectorRangeTreeMap &T, int parent, int child)
{
    auto diag = T.getRange(child);
    if (diag.second <= diag.first)
    {
        int pos = diag.first;
        diag = {pos, pos + 1}; // boundary placeholder
    }
    return DiagonalEdge{diag, parent, child};
}

static std::size_t partitionCacheCap()
{
    static std::size_t cap = 0;
    if (cap == 0)
    {
        if (const char *env = std::getenv("FLIPDIST_PARTITION_CACHE_CAP"))
        {
            long long v = std::atoll(env);
            if (v > 0)
                cap = static_cast<std::size_t>(v);
        }
        if (cap == 0)
            cap = 4000;
    }
    return cap;
}

static void maybeTrimPartitionCache()
{
    const std::size_t cap = partitionCacheCap();
    if (g_partitionSuccessMemo.size() > cap || g_partitionFailureMemo.size() > cap * 2)
    {
        g_partitionSuccessMemo.clear();
        g_partitionFailureMemo.clear();
    }
}

static std::string makePartitionEdgeCacheKey(const VectorRangeTreeMap &start,
                                             const VectorRangeTreeMap &target,
                                             const std::pair<int,int> &parent_range,
                                             const std::pair<int,int> &child_range)
{
    const auto &startSig = treeSignature(start);
    const auto &targetSig = treeSignature(target);
    std::string key;
    key.reserve(startSig.size() + targetSig.size() + 64);
    key.append(startSig);
    key.append("||");
    key.append(targetSig);
    key.append("|edge:");
    key.append(std::to_string(parent_range.first));
    key.push_back(',');
    key.append(std::to_string(parent_range.second));
    key.push_back('|');
    key.append(std::to_string(child_range.first));
    key.push_back(',');
    key.append(std::to_string(child_range.second));
    return key;
}

static std::string makePartitionRangeCacheKey(const VectorRangeTreeMap &start,
                                              const VectorRangeTreeMap &target,
                                              const std::pair<int,int> &diag)
{
    const auto &startSig = treeSignature(start);
    const auto &targetSig = treeSignature(target);
    std::string key;
    key.reserve(startSig.size() + targetSig.size() + 48);
    key.append(startSig);
    key.append("||");
    key.append(targetSig);
    key.append("|range:");
    key.append(std::to_string(diag.first));
    key.push_back(',');
    key.append(std::to_string(diag.second));
    return key;
}

static bool partitionAlongRangeCached(const VectorRangeTreeMap &start,
                                      const VectorRangeTreeMap &target,
                                      const std::pair<int,int> &diag,
                                      VectorRangeTreeMap &startInside,
                                      VectorRangeTreeMap &startOutside,
                                      VectorRangeTreeMap &targetInside,
                                      VectorRangeTreeMap &targetOutside);

// Pick a shared internal diagonal (by inorder range) that gives the most
// balanced split between its inside and outside parts. Returns false if none
static bool findBestCommonDiagonal(const VectorRangeTreeMap &start,
                                   const VectorRangeTreeMap &target,
                                   std::pair<int,int> &bestDiag)
{
    const int m = static_cast<int>(start.original_inorder.size());
    if (m <= 1 || static_cast<int>(target.original_inorder.size()) != m)
        return false;

    auto startIndex = buildEndpointIndex(start);
    auto targetIndex = buildEndpointIndex(target);
    if (startIndex.empty() || targetIndex.empty())
        return false;

    int bestScore = -1;
    std::pair<int,int> best{-1, -1};
    for (const auto &entry : startIndex)
    {
        const auto &diag = entry.first;
        if (!targetIndex.count(diag))
            continue;
        int width = diag.second - diag.first;
        if (width <= 1 || width >= m)
            continue;
        int inside = width;
        int outside = m - inside;
        if (inside <= 0 || outside <= 0)
            continue;
        int score = std::min(inside, outside);
        if (score > bestScore)
        {
            bestScore = score;
            best = diag;
        }
    }

    if (bestScore < 0)
        return false;
    bestDiag = best;
    return true;
}

// Repeatedly split both trees along shared diagonals, returning the leaf
// components after no further common diagonals remain
static bool decomposeAlongCommonDiagonals(
    const VectorRangeTreeMap &start,
    const VectorRangeTreeMap &target,
    std::vector<std::pair<VectorRangeTreeMap, VectorRangeTreeMap>> &components)
{
    components.clear();
    std::vector<std::pair<VectorRangeTreeMap, VectorRangeTreeMap>> stack;
    stack.emplace_back(start, target);
    bool split_any = false;

    while (!stack.empty())
    {
        auto current = std::move(stack.back());
        stack.pop_back();
        const auto &curStart = current.first;
        const auto &curTarget = current.second;

        std::pair<int,int> diag;
        if (!findBestCommonDiagonal(curStart, curTarget, diag))
        {
            components.push_back(std::move(current));
            continue;
        }

        VectorRangeTreeMap startInside;
        VectorRangeTreeMap startOutside;
        VectorRangeTreeMap targetInside;
        VectorRangeTreeMap targetOutside;
        if (!partitionAlongRangeCached(curStart, curTarget, diag,
                                       startInside, startOutside,
                                       targetInside, targetOutside))
        {
            components.push_back(std::move(current));
            continue;
        }

        if (startInside.original_nodes.empty() || startOutside.original_nodes.empty())
        {
            components.push_back(std::move(current));
            continue;
        }

        split_any = true;
        stack.emplace_back(std::move(startInside), std::move(targetInside));
        stack.emplace_back(std::move(startOutside), std::move(targetOutside));
    }

    return split_any && components.size() > 1;
}

static bool partitionAlongRangeCached(const VectorRangeTreeMap &start,
                                      const VectorRangeTreeMap &target,
                                      const std::pair<int,int> &diag,
                                      VectorRangeTreeMap &startInside,
                                      VectorRangeTreeMap &startOutside,
                                      VectorRangeTreeMap &targetInside,
                                      VectorRangeTreeMap &targetOutside)
{
    const std::string key = makePartitionRangeCacheKey(start, target, diag);
    if (g_partitionFailureMemo.count(key))
        return false;

    if (auto it = g_partitionSuccessMemo.find(key); it != g_partitionSuccessMemo.end())
    {
        const auto &entry = it->second;
        startInside = entry.startInside;
        startOutside = entry.startOutside;
        targetInside = entry.targetInside;
        targetOutside = entry.targetOutside;
        return true;
    }

    try
    {
        std::tie(startInside, startOutside) = VectorRangeTreeMap::partitionAlongRange(start, diag);
        std::tie(targetInside, targetOutside) = VectorRangeTreeMap::partitionAlongRange(target, diag);
    }
    catch (...)
    {
        g_partitionFailureMemo.insert(key);
        maybeTrimPartitionCache();
        return false;
    }

    if (startInside.original_nodes != targetInside.original_nodes ||
        startOutside.original_nodes != targetOutside.original_nodes)
    {
        g_partitionFailureMemo.insert(key);
        std::string keyRev = makePartitionRangeCacheKey(target, start, diag);
        g_partitionFailureMemo.insert(std::move(keyRev));
        maybeTrimPartitionCache();
        return false;
    }

    PartitionCacheEntry entry;
    entry.startInside = startInside;
    entry.startOutside = startOutside;
    entry.targetInside = targetInside;
    entry.targetOutside = targetOutside;
    entry.edgesInside = countInternalEdges(startInside);
    entry.edgesOutside = countInternalEdges(startOutside);
    g_partitionSuccessMemo.emplace(key, entry);

    std::string keyRev = makePartitionRangeCacheKey(target, start, diag);
    PartitionCacheEntry entryRev;
    entryRev.startInside = targetInside;
    entryRev.startOutside = targetOutside;
    entryRev.targetInside = startInside;
    entryRev.targetOutside = startOutside;
    entryRev.edgesInside = entry.edgesInside;
    entryRev.edgesOutside = entry.edgesOutside;
    g_partitionSuccessMemo.emplace(std::move(keyRev), std::move(entryRev));
    maybeTrimPartitionCache();
    return true;
}

enum class AdjSlot : int
{
    Parent = 0,
    Left = 1,
    Right = 2,
};

static DiagonalEdge makeBoundaryPlaceholderEdge(const VectorRangeTreeMap &T,
                                               int incidentNode,
                                               AdjSlot slot)
{
    auto r = T.getRange(incidentNode);
    int pos = r.first;
    if (slot == AdjSlot::Right)
        pos = r.second;
    // Use a degenerate interval (pos,pos) so this placeholder cannot collide
    // with any real diagonal interval (which always has L < R).
    const std::pair<int,int> diag{pos, pos};
    const int sentinelChild = VectorRangeTreeMap::NO_CHILD - 1 - static_cast<int>(slot);
    return DiagonalEdge{diag, incidentNode, sentinelChild};
}

static std::vector<DiagonalEdge> twoOtherIncidentEdgesOrBoundary(const VectorRangeTreeMap &T,
                                                                 int node,
                                                                 int neighbor)
{
    AdjSlot neighborSlot;
    if (T.getParent(node) == neighbor)
        neighborSlot = AdjSlot::Parent;
    else if (T.getLeftChild(node) == neighbor)
        neighborSlot = AdjSlot::Left;
    else if (T.getRightChild(node) == neighbor)
        neighborSlot = AdjSlot::Right;
    else
        return {};

    std::vector<DiagonalEdge> out;
    out.reserve(2);

    auto pushSlot = [&](AdjSlot slot) {
        switch (slot)
        {
            case AdjSlot::Parent:
            {
                int p = T.getParent(node);
                if (p != VectorRangeTreeMap::NO_PARENT && T.isOriginal(p))
                    out.push_back(makeDiagonalEdge(T, p, node));
                else
                    out.push_back(makeBoundaryPlaceholderEdge(T, node, slot));
                break;
            }
            case AdjSlot::Left:
            {
                int c = T.getLeftChild(node);
                if (c != VectorRangeTreeMap::NO_CHILD && T.isOriginal(c))
                    out.push_back(makeDiagonalEdge(T, node, c));
                else
                    out.push_back(makeBoundaryPlaceholderEdge(T, node, slot));
                break;
            }
            case AdjSlot::Right:
            {
                int c = T.getRightChild(node);
                if (c != VectorRangeTreeMap::NO_CHILD && T.isOriginal(c))
                    out.push_back(makeDiagonalEdge(T, node, c));
                else
                    out.push_back(makeBoundaryPlaceholderEdge(T, node, slot));
                break;
            }
        }
    };

    if (neighborSlot != AdjSlot::Parent)
        pushSlot(AdjSlot::Parent);
    if (neighborSlot != AdjSlot::Left)
        pushSlot(AdjSlot::Left);
    if (neighborSlot != AdjSlot::Right)
        pushSlot(AdjSlot::Right);

    if (out.size() != 2)
        return {};
    return out;
}

static bool findEdgeOrientation(const VectorRangeTreeMap &T,
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

static std::optional<std::pair<int,int>> orientDiagonalEdge(
        const VectorRangeTreeMap &T,
        const DiagonalEdge &opt,
        const std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> &diagMap)
{
    int pOrient, cOrient;
    // First try the stored orientation.
    if (findEdgeOrientation(T, opt.parent, opt.child, pOrient, cOrient))
        return std::make_pair(pOrient, cOrient);

    // Next try the current node that realises this diagonal (if present).
    auto it = diagMap.find(opt.diag);
    if (it == diagMap.end())
        return std::nullopt;

    int node = it->second;
    int parent = T.getParent(node);
    if (parent != VectorRangeTreeMap::NO_PARENT &&
        T.isOriginal(parent) &&
        findEdgeOrientation(T, parent, node, pOrient, cOrient))
    {
        return std::make_pair(pOrient, cOrient);
    }

    int left = T.getLeftChild(node);
    if (left != VectorRangeTreeMap::NO_CHILD &&
        T.isOriginal(left) &&
        findEdgeOrientation(T, node, left, pOrient, cOrient))
    {
        return std::make_pair(pOrient, cOrient);
    }

    int right = T.getRightChild(node);
    if (right != VectorRangeTreeMap::NO_CHILD &&
        T.isOriginal(right) &&
        findEdgeOrientation(T, node, right, pOrient, cOrient))
    {
        return std::make_pair(pOrient, cOrient);
    }

    // Final fallback: scan all nodes to find a child whose range matches the diagonal.
    for (int node : T.original_nodes)
    {
        auto r = T.getRange(node);
        if (r != opt.diag)
            continue;
        int parent = T.getParent(node);
        int p, c;
        if (parent != VectorRangeTreeMap::NO_PARENT && T.isOriginal(parent))
        {
            if (findEdgeOrientation(T, parent, node, p, c))
                return std::make_pair(p, c);
        }
        int left = T.getLeftChild(node);
        if (left != VectorRangeTreeMap::NO_CHILD && T.isOriginal(left))
        {
            if (findEdgeOrientation(T, node, left, p, c))
                return std::make_pair(p, c);
        }
        int right = T.getRightChild(node);
        if (right != VectorRangeTreeMap::NO_CHILD && T.isOriginal(right))
        {
            if (findEdgeOrientation(T, node, right, p, c))
                return std::make_pair(p, c);
        }
    }

    return std::nullopt;
}

static VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap &T)
{
    profile::ScopedTimer timer(&profile::g.ns_safeCopyTree);
    if (profile::enabled())
        ++profile::g.calls_safeCopyTree;

    try
    {
        // VectorRangeTreeMap is value-semantic; copying is safe and much cheaper
        // than rebuilding from traversals in tight recursion loops.
        return T;
    }
    catch (...)
    {
        // Return empty tree on failure.
    }
    return VectorRangeTreeMap{};
}

// Computes the diagonal (subtree range) that would be inserted by rotating the
// oriented edge (parent→child) in the given tree, without mutating the tree.
// In this BST-range model, only the ranges of the two rotated nodes change.
static std::optional<std::pair<int,int>> insertedDiagonalAfterRotate(
        const VectorRangeTreeMap &T,
        int parent,
        int child)
{
    if (!T.isOriginal(parent) || !T.isOriginal(child))
        return std::nullopt;

    auto itPos = T.position_in_inorder.find(parent);
    if (itPos == T.position_in_inorder.end())
        return std::nullopt;

    const int pos = itPos->second;

    // Right rotation at `parent` when `child` is the left child:
    // parent moves down and its new range is determined by (child.right, parent.right).
    if (T.getLeftChild(parent) == child)
    {
        int start = pos;
        int end   = pos + 1;

        int B = T.getRightChild(child);   // becomes parent.left
        int R = T.getRightChild(parent);  // unchanged
        if (B != VectorRangeTreeMap::NO_CHILD && T.isOriginal(B))
            start = T.getRange(B).first;
        if (R != VectorRangeTreeMap::NO_CHILD && T.isOriginal(R))
            end = T.getRange(R).second;
        return std::make_pair(start, end);
    }

    // Left rotation at `parent` when `child` is the right child:
    // parent moves down and its new range is determined by (parent.left, child.left).
    if (T.getRightChild(parent) == child)
    {
        int start = pos;
        int end   = pos + 1;

        int L = T.getLeftChild(parent);   // unchanged
        int B = T.getLeftChild(child);    // becomes parent.right
        if (L != VectorRangeTreeMap::NO_CHILD && T.isOriginal(L))
            start = T.getRange(L).first;
        if (B != VectorRangeTreeMap::NO_CHILD && T.isOriginal(B))
            end = T.getRange(B).second;
        return std::make_pair(start, end);
    }

    return std::nullopt;
}

// Free-edge search (no cache here; cache is in memo.cpp)
static std::vector<std::pair<int, int>> findFreeEdgesImpl(const VectorRangeTreeMap &T_init,
                                                         const VectorRangeTreeMap &T_final)
{
    profile::ScopedTimer timer(&profile::g.ns_findFreeEdges);
    if (profile::enabled())
        ++profile::g.calls_findFreeEdges;

    if (T_init.root == VectorRangeTreeMap::NO_CHILD ||
        T_final.root == VectorRangeTreeMap::NO_CHILD)
        return {};

    std::vector<std::pair<int, int>> candidates;
    auto targetDiagIndex = buildEndpointIndex(T_final);
    auto startDiagIndex  = buildEndpointIndex(T_init);

    std::unordered_set<std::pair<int,int>, PairHash, PairEq> initEdgesDirected;
    T_init.collectEdges(T_init.root, initEdgesDirected);

    for (const auto &edge : initEdgesDirected)
    {
        int parent = edge.first;
        int child  = edge.second;
        if (!T_init.isOriginal(parent) || !T_init.isOriginal(child))
            continue;

        // Identify the diagonal that would be inserted by rotating this edge.
        auto insertedDiagOpt = insertedDiagonalAfterRotate(T_init, parent, child);
        if (!insertedDiagOpt)
            continue;
        const auto insertedDiag = *insertedDiagOpt;
        if (insertedDiag.second - insertedDiag.first <= 1)
            continue;
        if (!targetDiagIndex.count(insertedDiag))
            continue;
        // Only consider rotations that actually insert a *new* target diagonal.
        if (startDiagIndex.count(insertedDiag))
            continue;

        candidates.push_back(edge); // keep original directed (parent->child)
    }

    return candidates;
}

static const std::vector<std::pair<int, int>> &findFreeEdgesCached(const VectorRangeTreeMap &T_init,
                                                                   const VectorRangeTreeMap &T_final)
{
    std::string key = makeFreeEdgeCacheKey(T_init, T_final);
    auto it = g_freeEdgeCache.find(key);
    if (it != g_freeEdgeCache.end())
        return it->second;

    auto edges = findFreeEdgesImpl(T_init, T_final);
    auto inserted = g_freeEdgeCache.emplace(std::move(key), std::move(edges));
    return inserted.first->second;
}

[[maybe_unused]] static std::string makeBranchMemoKey(std::size_t index, const std::unordered_set<int> &used_nodes)
{
    std::vector<int> nodes;
    nodes.reserve(used_nodes.size());
    for (int v : used_nodes)
        nodes.push_back(v);
    std::sort(nodes.begin(), nodes.end());
    std::string key = std::to_string(index) + ":";
    for (int v : nodes)
    {
        key += std::to_string(v);
        key.push_back(',');
    }
    return key;
}

struct BranchMemoKeyMask
{
    std::uint64_t used_mask = 0;
    std::uint64_t chosen_edges_lo = 0;
    std::uint64_t chosen_edges_hi = 0;
    std::uint32_t index = 0;
    bool operator==(const BranchMemoKeyMask &other) const noexcept
    {
        return used_mask == other.used_mask &&
               chosen_edges_lo == other.chosen_edges_lo &&
               chosen_edges_hi == other.chosen_edges_hi &&
               index == other.index;
    }
};

struct BranchMemoKeyMaskHash
{
    std::size_t operator()(const BranchMemoKeyMask &key) const noexcept
    {
        std::size_t seed = std::hash<std::uint64_t>{}(key.used_mask);
        seed ^= std::hash<std::uint64_t>{}(key.chosen_edges_lo) + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::uint64_t>{}(key.chosen_edges_hi) + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::uint32_t>{}(key.index) + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        return seed;
    }
};

static bool branchOnSPairsSet(const VectorRangeTreeMap &T_init,
                              const VectorRangeTreeMap &T_end,
                              int k,
                              const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S,
                              std::size_t index,
                              std::vector<std::pair<int,int>> &chosen,
                              std::unordered_set<int> &used_nodes,
                              std::unordered_set<std::string> &branchMemo,
                              const std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> &diagMap,
                              bool allow_independent_retry,
                              const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflictsPtr)
{
    if (index >= S.size())
    {
        if (chosen.empty())
            return false;
        return TreeDistI(T_init, T_end, k, chosen, allow_independent_retry, conflictsPtr);
    }

    if ((int)chosen.size() > k)
        return false;

    std::string memoKey;
    memoKey.reserve(64);
    memoKey.append(std::to_string(index));
    memoKey.push_back('|');
    memoKey.append(canonicalEdgeListKey(chosen));
    if (!branchMemo.insert(memoKey).second)
        return false;

    const auto &pair = S[index];
    const DiagonalEdge edges[2] = { pair.first, pair.second };

    // Try include-first or include-second if they’re independent.
    for (int i = 0; i < 2; ++i)
    {
        const auto &opt = edges[i];
        auto oriented = orientDiagonalEdge(T_init, opt, diagMap);
        if (!oriented)
        {
            // As a fallback, try the stored orientation directly if it is still an edge.
            if (hasParentChildEdge(T_init, opt.parent, opt.child))
            {
                oriented = std::make_pair(opt.parent, opt.child);
                debugPrint("branchOnSPairs: orientDiagonalEdge failed, using stored edge (" +
                           std::to_string(opt.parent) + "," + std::to_string(opt.child) + ")");
            }
            else
            {
                debugPrint("branchOnSPairs: unable to orient diagonal (" +
                           std::to_string(opt.diag.first) + "," + std::to_string(opt.diag.second) + "), skipping");
                continue;
            }
        }
        int parent = oriented->first;
        int child  = oriented->second;
        // Prefer branching only on edges whose diagonals are currently in conflict
        // (present in start but absent in target). This mirrors the Li–Xia
        // precondition C(Ts,Te)=0 after common-edge decomposition: common
        // diagonals never need to be flipped in an optimal sequence.
        if (conflictsPtr && conflictsPtr->count(opt.diag) == 0)
            continue;
        if (used_nodes.count(parent) || used_nodes.count(child))
            continue;

        used_nodes.insert(parent);
        used_nodes.insert(child);
        chosen.emplace_back(parent, child);
        if (branchOnSPairsSet(T_init, T_end, k, S, index + 1, chosen, used_nodes, branchMemo, diagMap, allow_independent_retry, conflictsPtr))
            return true;
        chosen.pop_back();
        used_nodes.erase(parent);
        used_nodes.erase(child);
    }

    // Branch include-none
    if (branchOnSPairsSet(T_init, T_end, k, S, index + 1, chosen, used_nodes, branchMemo, diagMap, allow_independent_retry, conflictsPtr))
        return true;

    return false;
}

static bool branchOnSPairsMask(const VectorRangeTreeMap &T_init,
                               const VectorRangeTreeMap &T_end,
                               int k,
                               const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S,
                               std::size_t index,
                               std::vector<std::pair<int,int>> &chosen,
                               std::uint64_t used_mask,
                               std::uint64_t chosen_edges_lo,
                               std::uint64_t chosen_edges_hi,
                               std::unordered_set<BranchMemoKeyMask, BranchMemoKeyMaskHash> &branchMemo,
                               const std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> &diagMap,
                               const std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> &edgeIdMap,
                               bool allow_independent_retry,
                               const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflictsPtr)
{
    if (index >= S.size())
    {
        if (chosen.empty())
            return false;
        return TreeDistI(T_init, T_end, k, chosen, allow_independent_retry, conflictsPtr);
    }

    if ((int)chosen.size() > k)
        return false;

    BranchMemoKeyMask memoKey{used_mask, chosen_edges_lo, chosen_edges_hi, static_cast<std::uint32_t>(index)};
    if (!branchMemo.insert(memoKey).second)
        return false;

    const auto &pair = S[index];
    const DiagonalEdge edges[2] = { pair.first, pair.second };

    // Try include-first or include-second if they’re independent.
    for (int i = 0; i < 2; ++i)
    {
        const auto &opt = edges[i];
        auto oriented = orientDiagonalEdge(T_init, opt, diagMap);
        if (!oriented)
        {
            if (hasParentChildEdge(T_init, opt.parent, opt.child))
            {
                oriented = std::make_pair(opt.parent, opt.child);
                debugPrint("branchOnSPairs: orientDiagonalEdge failed, using stored edge (" +
                           std::to_string(opt.parent) + "," + std::to_string(opt.child) + ")");
            }
            else
            {
                debugPrint("branchOnSPairs: unable to orient diagonal (" +
                           std::to_string(opt.diag.first) + "," + std::to_string(opt.diag.second) + "), skipping");
                continue;
            }
        }

        const int parent = oriented->first;
        const int child  = oriented->second;
        if (conflictsPtr && conflictsPtr->count(opt.diag) == 0)
            continue;
        if (parent < 0 || child < 0 || parent >= 64 || child >= 64)
            continue;
        const std::uint64_t bits = (1ull << parent) | (1ull << child);
        if (used_mask & bits)
            continue;

        auto itEdgeId = edgeIdMap.find({parent, child});
        if (itEdgeId == edgeIdMap.end())
            continue;
        int edgeId = itEdgeId->second;
        if (edgeId < 0 || edgeId >= 128)
            continue;
        std::uint64_t next_lo = chosen_edges_lo;
        std::uint64_t next_hi = chosen_edges_hi;
        if (edgeId < 64)
            next_lo |= (1ull << edgeId);
        else
            next_hi |= (1ull << (edgeId - 64));

        chosen.emplace_back(parent, child);
        if (branchOnSPairsMask(T_init, T_end, k, S, index + 1, chosen, used_mask | bits, next_lo, next_hi, branchMemo, diagMap, edgeIdMap, allow_independent_retry, conflictsPtr))
            return true;
        chosen.pop_back();
    }

    // Branch include-none
    if (branchOnSPairsMask(T_init, T_end, k, S, index + 1, chosen, used_mask, chosen_edges_lo, chosen_edges_hi, branchMemo, diagMap, edgeIdMap, allow_independent_retry, conflictsPtr))
        return true;

    return false;
}

static bool tryCommonEdgeDecomposition(const VectorRangeTreeMap &start,
                                       const VectorRangeTreeMap &target,
                                       int k)
{
    if (k <= 0)
        return false;

    // Decompose along any common polygon diagonal (identified by its inorder range),
    // not by matching parent/child labels. This mirrors the standard additivity
    // property for rotation distance along common internal diagonals.
    auto startIndex  = buildEndpointIndex(start);
    auto targetIndex = buildEndpointIndex(target);
    if (startIndex.empty() || targetIndex.empty())
        return false;

    const int m = static_cast<int>(start.original_inorder.size());
    if (m <= 1 || static_cast<int>(target.original_inorder.size()) != m)
        return false;

    struct Candidate
    {
        std::pair<int,int> diag;
        int balance = 0; // min(inside, outside)
        int width = 0;
    };

    std::vector<Candidate> candidates;
    candidates.reserve(std::min(startIndex.size(), targetIndex.size()));
    for (const auto &entry : startIndex)
    {
        const auto &diag = entry.first;
        if (!targetIndex.count(diag))
            continue;
        int width = diag.second - diag.first;
        int inside = width;
        int outside = m - inside;
        if (inside <= 0 || outside <= 0)
            continue;
        int balance = std::min(inside, outside);
        candidates.push_back(Candidate{diag, balance, width});
    }
    if (candidates.empty())
        return false;

    // Prefer balanced splits to shrink worst-case recursion; ties broken by
    // wider diagonals (often reduce conflicts more) and then lexicographically.
    std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) {
        if (a.balance != b.balance) return a.balance > b.balance;
        if (a.width != b.width) return a.width > b.width;
        return a.diag < b.diag;
    });

    // Try only a few best-balanced diagonals; this is a shortcut, not the main search.
    constexpr size_t COMMON_DIAG_TRY_CAP = 6;
    if (candidates.size() > COMMON_DIAG_TRY_CAP)
        candidates.resize(COMMON_DIAG_TRY_CAP);

    for (const auto &cand : candidates)
    {
        const auto &diag = cand.diag;
        try
        {
            VectorRangeTreeMap startInside;
            VectorRangeTreeMap startOutside;
            VectorRangeTreeMap targetInside;
            VectorRangeTreeMap targetOutside;
            if (!partitionAlongRangeCached(start, target, diag,
                                           startInside, startOutside,
                                           targetInside, targetOutside))
                continue;

            if (startInside.original_nodes.empty() || startOutside.original_nodes.empty())
                continue;

            int distInside = FlipDistMinK(startInside, targetInside, k);
            if (distInside < 0) continue;
            int remaining = k - distInside;
            if (remaining < 0) continue;
            int distOutside = FlipDistMinK(startOutside, targetOutside, remaining);
            if (distOutside < 0) continue;
            if (distInside + distOutside <= k)
                return true;
        }
        catch (...)
        {
            continue;
        }
    }

    return false;
}

static bool shouldForceLocalBFS(const VectorRangeTreeMap &start,
                                const VectorRangeTreeMap &target,
                                int k,
                                size_t pairCount)
{
    if (pairCount == 0 || k <= 0)
        return false;
    static std::unordered_set<std::string> visitedStates;

    std::string key;
    key.reserve(64);
    key += treeShapeSignature(start);
    key += "||";
    key += treeShapeSignature(target);
    key += "|k:";
    key += std::to_string(k);
    key += "|pairs:";
    key += std::to_string(pairCount);

    if (!visitedStates.insert(key).second)
    {
        if (visitedStates.size() > 4000)
            visitedStates.clear();
        return true;
    }
    if (visitedStates.size() > 4000)
        visitedStates.clear();
    return false;
}

static int minBudgetForTreeDistS(const VectorRangeTreeMap &start,
                                 const VectorRangeTreeMap &target,
                                 const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S_side,
                                 int maxBudget,
                                 bool allow_independent_retry)
{
    int lb = lowerBoundEdgeDifference(start, target);
    if (lb > maxBudget)
        return -1;

    // Consult TreeDistS's monotone bounds cache first.
    const auto &startSig = treeShapeSignature(start);
    const auto &targetSig = treeShapeSignature(target);
    std::string baseKey;
    baseKey.reserve(startSig.size() + targetSig.size() + 2);
    baseKey.append(startSig);
    baseKey.append("||");
    baseKey.append(targetSig);
    if (auto itBounds = g_treeDistSBounds.find(baseKey); itBounds != g_treeDistSBounds.end())
    {
        const auto &bounds = itBounds->second;
        if (bounds.min_true != 0x3fffffff && bounds.min_true <= maxBudget)
            return bounds.min_true;
        if (bounds.max_false >= maxBudget)
            return -1;
    }

    // If the lower bound itself works, we're done.
    if (TreeDistS(start, target, lb, S_side, allow_independent_retry, nullptr))
        return lb;

    // Exponential search to find some satisfiable budget, then binary search
    // for the minimum. This is significantly cheaper than linear scanning
    // when k is large and TreeDistS is invoked repeatedly by partitioning.
    int lo = lb;
    int step = 1;
    int hi = -1;
    while (lo + step <= maxBudget)
    {
        int probe = lo + step;
        if (TreeDistS(start, target, probe, S_side, allow_independent_retry, nullptr))
        {
            hi = probe;
            break;
        }
        lo = probe;
        step *= 2;
    }
    if (hi < 0)
    {
        if (!TreeDistS(start, target, maxBudget, S_side, allow_independent_retry, nullptr))
            return -1;
        hi = maxBudget;
    }

    int left = lo + 1;
    int right = hi;
    while (left < right)
    {
        int mid = left + (right - left) / 2;
        if (TreeDistS(start, target, mid, S_side, allow_independent_retry, nullptr))
            right = mid;
        else
            left = mid + 1;
    }
    return left;
}

// Explore independent sets recursively.
static bool enumerateIndependentSetsRecursive(
        const VectorRangeTreeMap &start,
        const VectorRangeTreeMap &target,
        int k,
        const std::vector<std::pair<int,int>> &edges,
        size_t index,
        std::vector<std::pair<int,int>> &current,
        std::vector<char> &used,
        bool allow_independent_retry)
{
    if ((int)current.size() > k)
        return false;

    if (index == edges.size())
    {
        if (current.empty())
            return false;
        return TreeDistI(start, target, k, current, allow_independent_retry, nullptr);
    }

    const auto &edge = edges[index];

    // Try including the current edge first so we prioritise larger independent sets.
    if (edge.first >= 0 && edge.first < (int)used.size() &&
        edge.second >= 0 && edge.second < (int)used.size() &&
        !used[edge.first] && !used[edge.second])
    {
        used[edge.first]  = 1;
        used[edge.second] = 1;
        current.push_back(edge);
        if (enumerateIndependentSetsRecursive(start, target, k, edges, index + 1, current, used, allow_independent_retry))
            return true;
        current.pop_back();
        used[edge.first]  = 0;
        used[edge.second] = 0;
    }

    // Exclude current edge
    return enumerateIndependentSetsRecursive(start, target, k, edges, index + 1, current, used, allow_independent_retry);
}

static bool exploreIndependentSets(const VectorRangeTreeMap &start,
                                   const VectorRangeTreeMap &target,
                                   int k,
                                   const std::vector<std::pair<int,int>> &edges,
                                   bool allow_independent_retry)
{
    std::vector<std::pair<int,int>> current;
    current.reserve(edges.size());
    int usedSize = std::max(start.max_node_value, target.max_node_value) + 2;
    size_t allocSize = usedSize > 0 ? static_cast<size_t>(usedSize) : 0;
    std::vector<char> used(allocSize, 0);

    return enumerateIndependentSetsRecursive(start, target, k, edges, 0, current, used, allow_independent_retry);
}

// Ear contraction (optional)
static bool contractSharedEar(VectorRangeTreeMap &start,
                              VectorRangeTreeMap &target,
                              const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflicts)
{
    if (start.root == VectorRangeTreeMap::NO_CHILD ||
        target.root == VectorRangeTreeMap::NO_CHILD)
        return false;

    auto startIndex  = buildEndpointIndex(start);
    auto targetIndex = buildEndpointIndex(target);

    auto skipDiag = [&](const std::pair<int,int> &diag) {
        return conflicts && conflicts->count(diag);
    };

    for (const auto &entry : startIndex)
    {
        const auto &key = entry.first; // (L,R)
        auto itTarget = targetIndex.find(key);
        if (itTarget == targetIndex.end()) continue;
        if (skipDiag(key))
            continue;
        if (conflicts && conflicts->count(key))
            continue;

        int nodeStart  = entry.second;
        int nodeTarget = itTarget->second;

        auto parentRangeStart  = start.getRange(nodeStart);
        auto parentRangeTarget = target.getRange(nodeTarget);

        int children[2] = { start.getLeftChild(nodeStart), start.getRightChild(nodeStart) };

        for (int child : children)
        {
            if (child == VectorRangeTreeMap::NO_CHILD) continue;
            if (!start.isOriginal(child)) continue;

            auto childRange = start.getRange(child);
            if (childRange.first >= childRange.second) continue;

            auto startParts  = VectorRangeTreeMap::partitionAlongEdge(start,  parentRangeStart,  childRange);
            auto targetParts = VectorRangeTreeMap::partitionAlongEdge(target, parentRangeTarget, childRange);

            const auto &startChildTree   = startParts.first;
            const auto &startComplement  = startParts.second;
            const auto &targetChildTree  = targetParts.first;
            const auto &targetComplement = targetParts.second;

            if (startChildTree.original_nodes == targetChildTree.original_nodes &&
                TreesEqual(startChildTree, targetChildTree))
            {
                start  = startComplement;
                target = targetComplement;
                return true;
            }

            if (startComplement.original_nodes == targetComplement.original_nodes &&
                TreesEqual(startComplement, targetComplement))
            {
                start  = startChildTree;
                target = targetChildTree;
                return true;
            }
        }
    }

    for (const auto &entry : startIndex)
    {
        const auto &diag = entry.first;
        auto itTarget = targetIndex.find(diag);
        if (itTarget == targetIndex.end()) continue;
        if (skipDiag(diag))
            continue;

        try
        {
            VectorRangeTreeMap startInside;
            VectorRangeTreeMap startOutside;
            VectorRangeTreeMap targetInside;
            VectorRangeTreeMap targetOutside;
            if (!partitionAlongRangeCached(start, target, diag,
                                           startInside, startOutside,
                                           targetInside, targetOutside))
                continue;

            if (startInside.original_nodes == targetInside.original_nodes &&
                TreesEqual(startInside, targetInside))
            {
                start  = startOutside;
                target = targetOutside;
                return true;
            }

            if (startOutside.original_nodes == targetOutside.original_nodes &&
                TreesEqual(startOutside, targetOutside))
            {
                start  = startInside;
                target = targetInside;
                return true;
            }
        }
        catch (...)
        {
            continue;
        }
    }

    return false;
}

} // namespace

bool hasParentChildEdge(const VectorRangeTreeMap &T, int parent, int child)
{
    try
    {
        if (!T.isOriginal(parent) || !T.isOriginal(child))
            return false;
        int left = T.getLeftChild(parent);
        int right = T.getRightChild(parent);
        return (left == child) || (right == child);
    }
    catch (...)
    {
        return false;
    }
}

bool FlipDistTree(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k)
{
    profile::ScopedTimer timer(&profile::g.ns_flipDistTree);
    if (profile::enabled())
        ++profile::g.calls_flipDistTree;
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    if (k < 0)
        return false;

    std::string startKey = treeShapeSignature(T_init);
    std::string targetKey = treeShapeSignature(T_final);
    FlipMemoKey memoKey{startKey, targetKey, k};
    FlipMemoKey memoKeyRev{targetKey, startKey, k};
    std::string baseKey = startKey;
    baseKey.append("||");
    baseKey.append(targetKey);
    std::string baseKeyRev = targetKey;
    baseKeyRev.append("||");
    baseKeyRev.append(startKey);
    if (auto it = g_flipDistMemo.find(memoKey); it != g_flipDistMemo.end())
    {
        debugPrint(std::string("FlipDistTree memo hit k=") + std::to_string(k));
        if (profile::enabled())
            ++profile::g.memo_hits_flipDistTree;
        return it->second;
    }
    if (auto it = g_flipDistMemo.find(memoKeyRev); it != g_flipDistMemo.end())
    {
        debugPrint(std::string("FlipDistTree memo hit (reversed) k=") + std::to_string(k));
        if (profile::enabled())
            ++profile::g.memo_hits_flipDistTree;
        return it->second;
    }

    auto itBounds = g_flipDistBounds.find(baseKey);
    if (itBounds != g_flipDistBounds.end())
    {
        const auto &bounds = itBounds->second;
        if (bounds.max_false >= 0 && k <= bounds.max_false)
            return false;
        if (bounds.min_true != 0x3fffffff && k >= bounds.min_true)
            return true;
    }
    auto itBoundsRev = g_flipDistBounds.find(baseKeyRev);
    if (itBoundsRev != g_flipDistBounds.end())
    {
        const auto &bounds = itBoundsRev->second;
        if (bounds.max_false >= 0 && k <= bounds.max_false)
            return false;
        if (bounds.min_true != 0x3fffffff && k >= bounds.min_true)
            return true;
    }

    int lb_conflicts = lowerBoundEdgeDifference(T_init, T_final);
    if (earContractionEnabledFor(lb_conflicts))
    {
        VectorRangeTreeMap reducedStart = T_init;
        VectorRangeTreeMap reducedTarget = T_final;
        bool contracted = false;
        while (contractSharedEar(reducedStart, reducedTarget, nullptr))
        {
            contracted = true;
        }
        if (contracted)
        {
            return FlipDistTree(reducedStart, reducedTarget, k);
        }
    }

    auto memoReturn = [&](bool result) {
        g_flipDistMemo[memoKey] = result;
        g_flipDistMemo[memoKeyRev] = result;
        auto &bounds = g_flipDistBounds[baseKey];
        auto &boundsRev = g_flipDistBounds[baseKeyRev];
        if (result)
        {
            bounds.min_true = std::min(bounds.min_true, k);
            boundsRev.min_true = std::min(boundsRev.min_true, k);
        }
        else
        {
            bounds.max_false = std::max(bounds.max_false, k);
            boundsRev.max_false = std::max(boundsRev.max_false, k);
        }
        return result;
    };

    bool solved = false;
    auto tryRecord = [&](bool v) {
        if (v) solved = true;
    };

    auto conflicts = collectConflictingEdges(T_init, T_final);
    if (conflicts.empty())
    {
        debugPrint("FlipDistTree: conflicts empty");
        // No conflicts: trees share all diagonals.
        bool result = countInternalEdges(T_final) == 0;
        return memoReturn(result);
    }
    else
    {
        debugPrint("FlipDistTree: conflicts size=" + std::to_string(conflicts.size()));
    }

    // PSEUDOCODE STEP 0: lower bound check
    if (lb_conflicts > k)
    {
        return memoReturn(false);
    }

    // If there is already a free edge (a rotation that inserts a target
    // diagonal), TreeDistS should process it immediately (Li–Xia step 1)
    // rather than spending time enumerating initial independent sets.
    if (k > 0)
    {
        const auto &freeEdges = findFreeEdgesCached(T_init, T_final);
        if (!freeEdges.empty())
        {
            bool result = TreeDistS(T_init, T_final, k, {}, /*allow_independent_retry=*/true, nullptr);
            return memoReturn(result);
        }
    }

    // PSEUDOCODE STEP 1: Enumerate independent internal edge sets I in T_init.
    // (In this representation, internal edges are the parent→child edges.)
    // This keeps the solver complete with respect to Li–Xia’s top-level
    // branching while still being tractable for n<=25 (the number of
    // independent sets over a tree’s edges is manageable).
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> conflictEdgeSet;
    conflictEdgeSet.reserve(conflicts.size());
    for (const auto &c : conflicts)
        conflictEdgeSet.insert({c.parent, c.child});

    std::vector<std::pair<int,int>> candidateEdges = getInternalEdges(T_init);

    std::sort(candidateEdges.begin(), candidateEdges.end(), [&](const auto &a, const auto &b) {
        bool ca = conflictEdgeSet.count(a) > 0;
        bool cb = conflictEdgeSet.count(b) > 0;
        if (ca != cb) return ca > cb;
        auto da = T_init.diagonalEndpoints(a.second);
        auto db = T_init.diagonalEndpoints(b.second);
        int wa = da.second - da.first;
        int wb = db.second - db.first;
        if (wa != wb) return wa > wb;
        return a < b;
    });

    // First try independent sets restricted to conflicting edges only.
    std::vector<std::pair<int,int>> conflictEdges;
    conflictEdges.reserve(candidateEdges.size());
    for (const auto &edge : candidateEdges)
    {
        if (conflictEdgeSet.count(edge))
            conflictEdges.push_back(edge);
    }
    if (!conflictEdges.empty())
    {
        tryRecord(exploreIndependentSets(T_init, T_final, k, conflictEdges, true));
        if (solved)
            return memoReturn(true);
    }

    // Next, try independent sets over the conflict neighborhood (incident edges).
    std::vector<std::pair<int,int>> neighborEdges;
    if (!conflictEdges.empty())
    {
        std::unordered_set<std::pair<int,int>, PairHash, PairEq> neighborSet;
        neighborSet.reserve(conflictEdges.size() * 3);
        for (const auto &edge : conflictEdges)
        {
            neighborSet.insert(edge);
            auto incParent = getIncidentEdges(T_init, edge.first);
            auto incChild = getIncidentEdges(T_init, edge.second);
            for (const auto &inc : incParent)
                neighborSet.insert(inc);
            for (const auto &inc : incChild)
                neighborSet.insert(inc);
        }
        neighborEdges.assign(neighborSet.begin(), neighborSet.end());
        std::sort(neighborEdges.begin(), neighborEdges.end(), [&](const auto &a, const auto &b) {
            bool ca = conflictEdgeSet.count(a) > 0;
            bool cb = conflictEdgeSet.count(b) > 0;
            if (ca != cb) return ca > cb;
            auto da = T_init.diagonalEndpoints(a.second);
            auto db = T_init.diagonalEndpoints(b.second);
            int wa = da.second - da.first;
            int wb = db.second - db.first;
            if (wa != wb) return wa > wb;
            return a < b;
        });

        if (neighborEdges.size() > conflictEdges.size() &&
            neighborEdges.size() < candidateEdges.size())
        {
            tryRecord(exploreIndependentSets(T_init, T_final, k, neighborEdges, true));
            if (solved)
                return memoReturn(true);
        }
    }

    bool allow_full_enum = true;
    static int full_conflict_cap = -2;
    if (full_conflict_cap == -2)
    {
        if (const char *env = std::getenv("FLIPDIST_I_FULL_CONFLICT_CAP"))
            full_conflict_cap = std::max(0, std::atoi(env));
        else
            full_conflict_cap = -1;
    }
    if (full_conflict_cap >= 0 && static_cast<int>(conflictEdges.size()) > full_conflict_cap)
    {
        allow_full_enum = false;
        debugPrint("FlipDistTree: skipping full independent-set enumeration (conflictEdges=" +
                   std::to_string(conflictEdges.size()) +
                   ", cap=" + std::to_string(full_conflict_cap) + ")");
    }

    // Fall back to full independent-set enumeration for completeness.
    if (allow_full_enum)
        tryRecord(exploreIndependentSets(T_init, T_final, k, candidateEdges, true));
    if (solved)
    {
        return memoReturn(true);
    }

    // Global fallback: if budget is small, run a capped BFS to retain completeness.
    if (!solved && k <= 20 && shouldUseInternalBFS(T_init, k))
    {
        int bfsCap = k;
        int bfsDist = profiledMinRotationsBFS(T_init, T_final, bfsCap);
        if (bfsDist >= 0 && bfsDist <= k)
        {
            debugPrint("FlipDistTree: global BFS succeeded with dist=" + std::to_string(bfsDist));
            return memoReturn(true);
        }
    }

    // No solution found within k
    if (!solved)
        debugPrint("No solution found, returning false");
    return memoReturn(solved);
}

int FlipDistMinK(const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2, int k_max)
{
    if (TreesEqual(T1, T2))
        return 0;
    int lb0 = lowerBoundEdgeDifference(T1, T2);
    if (lb0 > k_max)
        return -1;

    // If the trees share internal diagonals, rotation distance is additive
    // across those splits. Fully decompose along common diagonals first to
    // shrink the problem before invoking Li–Xia branching.
    std::vector<std::pair<VectorRangeTreeMap, VectorRangeTreeMap>> components;
    if (decomposeAlongCommonDiagonals(T1, T2, components))
    {
        std::sort(components.begin(), components.end(), [](const auto &a, const auto &b) {
            return countInternalEdges(a.first) < countInternalEdges(b.first);
        });

        int total = 0;
        for (const auto &comp : components)
        {
            int remaining = k_max - total;
            if (remaining < 0)
                return -1;
            int dist = FlipDistMinK(comp.first, comp.second, remaining);
            if (dist < 0)
                return -1;
            total += dist;
            if (total > k_max)
                return -1;
        }
        return total;
    }

    // Outer k-search: exponential search to find a satisfiable k, then binary
    // search for the minimum. This avoids spending O(dist-lb) full runs proving
    // failure when the true distance is far above the lower bound.
    //
    // We still start at the lower bound (never below budget), preserving the
    // Li–Xia FPT behaviour inside each FlipDistTree call.
    int lo = std::max(lb0, 0);
    if (lo == 0 && FlipDistTree(T1, T2, 0))
        return 0;
    if (lo > 0 && FlipDistTree(T1, T2, lo))
        return lo;

    int step = 1;
    int hi = -1;
    int last_false = lo;
    for (;;)
    {
        int probe = lo + step;
        if (probe > k_max)
            break;
        if (FlipDistTree(T1, T2, probe))
        {
            hi = probe;
            break;
        }
        last_false = probe;
        step *= 2;
        lo = probe;
    }
    if (hi < 0)
        return -1;

    int left = last_false + 1;
    int right = hi;
    while (left < right)
    {
        int mid = left + (right - left) / 2;
        if (FlipDistTree(T1, T2, mid))
            right = mid;
        else
            left = mid + 1;
    }
    return left;
}

// TreeDistI: Handles independent edge set I
bool TreeDistI(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k,
               const std::vector<std::pair<int, int>> &I,
               bool allow_independent_retry,
               const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflictSet)
{
    profile::ScopedTimer timer(&profile::g.ns_treeDistI);
    if (profile::enabled())
        ++profile::g.calls_treeDistI;
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    // TreeDistI is frequently called with many distinct I-sets; memoizing large
    // independent sets tends to be low-hit and expensive due to key construction.
    constexpr size_t TREE_DIST_I_MEMO_MAX_I = 6;

    std::string memoKey;
    if (!I.empty() && I.size() <= TREE_DIST_I_MEMO_MAX_I)
    {
        memoKey.reserve(200);
        // Use shape signatures so identical induced subproblems (up to inorder relabeling)
        // share memo entries. TreeDistI gets called extremely often inside TreeDistS,
        // and label-specific keys miss most reuse.
        memoKey += treeShapeSignature(T_init);
        memoKey += "||";
        memoKey += treeShapeSignature(T_final);
        memoKey += "|k:";
        memoKey += std::to_string(k);
        memoKey += "|I:";
        // Encode the chosen oriented edges in the canonical inorder-rank space.
        std::vector<std::pair<int,int>> rankEdges;
        rankEdges.reserve(I.size());
        for (const auto &edge : I)
        {
            int parent = edge.first;
            int child  = edge.second;
            auto itP = T_init.position_in_inorder.find(parent);
            auto itC = T_init.position_in_inorder.find(child);
            if (itP == T_init.position_in_inorder.end() || itC == T_init.position_in_inorder.end())
            {
                rankEdges.clear();
                break;
            }
            rankEdges.emplace_back(itP->second + 1, itC->second + 1);
        }
        if (rankEdges.size() != I.size())
        {
            memoKey.clear();
        }
        else
        {
            memoKey += canonicalEdgeListKey(rankEdges);
        }
        if (!memoKey.empty())
        {
            auto it = g_treeDistIMemo.find(memoKey);
            if (it != g_treeDistIMemo.end())
            {
                if (profile::enabled())
                    ++profile::g.memo_hits_treeDistI;
                debugPrint("TreeDistI: memo hit");
                return it->second;
            }
        }
    }

    int remaining_budget = k - (int)I.size(); // k - |I| from pseudocode

    if (remaining_budget < 0)
    {
        debugPrint("TreeDistI: Not enough budget for |I| rotations");
        return false;
    }

    int lb_conflicts = (conflictSet) ? (int)conflictSet->size() : lowerBoundEdgeDifference(T_init, T_final);
    debugPrint("TreeDistI: lower bound=" + std::to_string(lb_conflicts) +
               " remaining_budget=" + std::to_string(remaining_budget) +
               " |I|=" + std::to_string(I.size()));
    // Be permissive: do not prune here; later recursion and k bounds will enforce budget.

    // Special handling when budget exactly equals |I|
    if (remaining_budget == 0)
    {
        VectorRangeTreeMap T_bar = safeCopyTree(T_init);
        for (const auto &edge : I)
        {
            int parent = edge.first, child = edge.second;
            if (!hasParentChildEdge(T_bar, parent, child))
            {
                int pOrient, cOrient;
                if (findEdgeOrientation(T_bar, parent, child, pOrient, cOrient))
                {
                    parent = pOrient;
                    child  = cOrient;
                }
                else
                {
                    debugPrint("TreeDistI: Invalid edge in budget==0 path; rejecting branch");
                    if (!memoKey.empty())
                        g_treeDistIMemo[memoKey] = false;
                    return false;
                }
            }

            if (T_bar.getLeftChild(parent) == child)
                T_bar.rotateRight(parent);
            else if (T_bar.getRightChild(parent) == child)
                T_bar.rotateLeft(parent);
            else
            {
                debugPrint("TreeDistI: Edge not adjacent in budget==0 path; rejecting branch");
                if (!memoKey.empty())
                    g_treeDistIMemo[memoKey] = false;
                return false;
            }
        }

        bool equal = TreesEqual(T_bar, T_final);
    if (!memoKey.empty())
        g_treeDistIMemo[memoKey] = equal;
        return equal;
    }

    // PSEUDOCODE STEP 2: Rotate edges in I and build S
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S;

    VectorRangeTreeMap T_bar = safeCopyTree(T_init);
    if (T_bar.original_nodes.empty())
    {
        debugPrint("TreeDistI: Failed to copy tree");
        return false;
    }

    for (const auto &edge : I)
    {
        int parent = edge.first;
        int child = edge.second;

        debugPrint("TreeDistI: Processing edge (" + std::to_string(parent) + "," + std::to_string(child) + ")");

        if (!hasParentChildEdge(T_bar, parent, child))
        {
            int pOrient, cOrient;
            if (findEdgeOrientation(T_bar, parent, child, pOrient, cOrient))
            {
                parent = pOrient;
                child  = cOrient;
            }
            else
            {
                debugPrint("TreeDistI: Invalid edge; rejecting branch");
                if (!memoKey.empty())
                    g_treeDistIMemo[memoKey] = false;
                return false;
            }
        }

        int u, v; // nodes joined by new edge

        if (T_bar.getLeftChild(parent) == child)
        {
            u = child;
            v = parent;
            T_bar.rotateRight(parent);
        }
        else if (T_bar.getRightChild(parent) == child)
        {
            u = child;
            v = parent;
            T_bar.rotateLeft(parent);
        }
        else
        {
            debugPrint("TreeDistI: Edge not adjacent; rejecting branch");
            if (!memoKey.empty())
                g_treeDistIMemo[memoKey] = false;
            return false;
        }

        auto u_pair = twoOtherIncidentEdgesOrBoundary(T_bar, u, v);
        if (u_pair.size() == 2)
        {
            if (appendValidPair(S, u_pair[0], u_pair[1]))
                debugPrint("TreeDistI: Added edge pair for u");
        }

        auto v_pair = twoOtherIncidentEdgesOrBoundary(T_bar, v, u);
        if (v_pair.size() == 2)
        {
            if (appendValidPair(S, v_pair[0], v_pair[1]))
                debugPrint("TreeDistI: Added edge pair for v");
        }
    }

    deduplicateDiagonalPairs(S, "TreeDistI_S");

    // Optional pruning: require TreeDistI rotations to not "worsen" (or strictly
    // improve) the conflict-based lower bound before recursing into TreeDistS.
    // This is disabled by default for completeness; enable via
    // `FLIPDIST_CONFLICT_PRUNE`:
    //   0 = off (default)
    //   1 = reject if bound increases
    //   2 = reject if bound does not decrease
    static int prune_mode = -1;
    if (prune_mode == -1)
    {
        if (const char *env = std::getenv("FLIPDIST_CONFLICT_PRUNE"))
            prune_mode = std::atoi(env);
        else
            prune_mode = 0;
    }
    if (prune_mode > 0)
    {
        int bound_before = lowerBoundEdgeDifference(T_init, T_final);
        int bound_after  = lowerBoundEdgeDifference(T_bar, T_final);
        if ((prune_mode == 1 && bound_after > bound_before) ||
            (prune_mode >= 2 && bound_after >= bound_before))
        {
            debugPrint("TreeDistI: conflict-prune bound " + std::to_string(bound_before) +
                       " -> " + std::to_string(bound_after) + ", rejecting branch");
            if (!memoKey.empty())
                g_treeDistIMemo[memoKey] = false;
            return false;
        }
    }

    if (DEBUG)
    {
        std::string msg = "TreeDistI: S pairs:";
        if (S.empty())
        {
            msg += " <empty>";
        }
        else
        {
            for (const auto &pair : S)
            {
                msg += " (diag(" + std::to_string(pair.first.diag.first) + "," + std::to_string(pair.first.diag.second) +
                       ") edge(" + std::to_string(pair.first.parent) + "," + std::to_string(pair.first.child) + ");" +
                       " diag(" + std::to_string(pair.second.diag.first) + "," + std::to_string(pair.second.diag.second) +
                       ") edge(" + std::to_string(pair.second.parent) + "," + std::to_string(pair.second.child) + "))";
            }
        }
        debugPrint(msg);
    }

    // Safe pruning: if the lower bound on the rotated tree already exceeds
    // the remaining budget, this branch cannot succeed.
    int lb_after = lowerBoundEdgeDifference(T_bar, T_final);
    if (lb_after > remaining_budget)
    {
        debugPrint("TreeDistI: lb_after=" + std::to_string(lb_after) +
                   " exceeds remaining_budget=" + std::to_string(remaining_budget));
        if (!memoKey.empty())
            g_treeDistIMemo[memoKey] = false;
        return false;
    }

    // Previously we pruned on non-improving conflict counts; for completeness, skip that here.

    bool result = TreeDistS(T_bar, T_final, k - (int)I.size(), S, allow_independent_retry, conflictSet);
    if (!memoKey.empty())
        g_treeDistIMemo[memoKey] = result;
    return result;
}

// TreeDistS: Handles S-branching and partitioning
bool TreeDistS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
               const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S,
               bool allow_independent_retry,
               const std::unordered_set<std::pair<int,int>, PairHash, PairEq> *conflictsPtr)
{
    profile::ScopedTimer timer(&profile::g.ns_treeDistS);
    if (profile::enabled())
        ++profile::g.calls_treeDistS;
    profile::heartbeat();

    constexpr bool USE_S_MEMO = true;
    static int conflict_dump_budget = 5;
    static int partner_dump_budget = 8;
    static int trace_budget = 8; // small targeted trace to understand stalls
    debugPrint("TreeDistS: start k=" + std::to_string(k) + " S size=" + std::to_string(S.size()) +
               " conflicts=" + std::to_string(conflictsPtr ? conflictsPtr->size() : 0));
    if (DEBUG && trace_budget > 0)
    {
        --trace_budget;
        std::ostringstream oss;
        oss << "[TRACE] TreeDistS k=" << k
            << " conflicts=" << (conflictsPtr ? conflictsPtr->size() : lowerBoundEdgeDifference(T_init, T_end))
            << " S_size=" << S.size();
        debugPrint(oss.str());
    }

    std::string memoKey;
    std::string baseKey;
    auto memoReturn = [&](bool result) -> bool {
        if (USE_S_MEMO && !memoKey.empty())
            g_treeDistSMemo[memoKey] = result;
        if (USE_S_MEMO && !baseKey.empty())
        {
            auto &bounds = g_treeDistSBounds[baseKey];
            if (result)
                bounds.min_true = std::min(bounds.min_true, k);
            else
                bounds.max_false = std::max(bounds.max_false, k);
            static std::size_t boundsCap = 0;
            if (boundsCap == 0)
            {
                if (const char *env = std::getenv("FLIPDIST_BOUNDS_CAP"))
                {
                    long long v = std::atoll(env);
                    if (v > 0)
                        boundsCap = static_cast<std::size_t>(v);
                }
                if (boundsCap == 0)
                    boundsCap = 2'000'000;
            }
            if (g_treeDistSBounds.size() > boundsCap)
                g_treeDistSBounds.clear();
        }
        return result;
    };

    if (USE_S_MEMO)
    {
        const auto &startSig = treeShapeSignature(T_init);
        const auto &endSig = treeShapeSignature(T_end);
        baseKey.clear();
        baseKey.reserve(startSig.size() + endSig.size() + 2);
        baseKey.append(startSig);
        baseKey.append("||");
        baseKey.append(endSig);

        auto checkBounds = [&](const std::string &key) -> std::optional<bool> {
            auto it = g_treeDistSBounds.find(key);
            if (it == g_treeDistSBounds.end())
                return std::nullopt;
            const auto &bounds = it->second;
            if (bounds.max_false >= 0 && k <= bounds.max_false)
                return false;
            if (bounds.min_true != 0x3fffffff && k >= bounds.min_true)
                return true;
            return std::nullopt;
        };

        if (auto bound = checkBounds(baseKey); bound.has_value())
            return *bound;

        std::string baseKeyRev = endSig;
        baseKeyRev.append("||");
        baseKeyRev.append(startSig);
        if (auto bound = checkBounds(baseKeyRev); bound.has_value())
            return *bound;

        memoKey = baseKey;
        memoKey.reserve(memoKey.size() + 16);
        memoKey.append("|k:");
        memoKey.append(std::to_string(k));
        if (auto itMemo = g_treeDistSMemo.find(memoKey); itMemo != g_treeDistSMemo.end())
        {
            if (profile::enabled())
                ++profile::g.memo_hits_treeDistS;
            debugPrint("TreeDistS: memo hit");
            return itMemo->second;
        }
    }

    if (k < 0)
        return memoReturn(false);

    if (TreesEqual(T_init, T_end))
        return memoReturn(true);

    int lb_conflicts = lowerBoundEdgeDifference(T_init, T_end);
    debugPrint("TreeDistS: lb_conflicts=" + std::to_string(lb_conflicts) + " k=" + std::to_string(k));
    if (lb_conflicts > k)
        return memoReturn(false);
    if (profile::enabled())
        profile::g.max_conflicts_seen = std::max<std::uint64_t>(profile::g.max_conflicts_seen, static_cast<std::uint64_t>(lb_conflicts));

    if (S.empty() && earContractionEnabledFor(lb_conflicts))
    {
        VectorRangeTreeMap reducedStart = T_init;
        VectorRangeTreeMap reducedTarget = T_end;
        bool contracted = false;
        while (contractSharedEar(reducedStart, reducedTarget, nullptr))
            contracted = true;
        if (contracted)
            return memoReturn(TreeDistS(reducedStart, reducedTarget, k, S, allow_independent_retry, nullptr));
    }

    // Try common edge decomposition shortcut (only when S is empty).
    if (S.empty() && tryCommonEdgeDecomposition(T_init, T_end, k))
        return memoReturn(true);

    bool attemptedFreeEdge = false;
    if (k > 0)
    {
        attemptedFreeEdge = true;
    const auto &freeEdges = findFreeEdgesCached(T_init, T_end);
    if (profile::enabled())
        profile::g.total_freeEdgeCandidates += static_cast<std::uint64_t>(freeEdges.size());
    debugPrint("TreeDistS: free edge candidates = " + std::to_string(freeEdges.size()));
    if (!freeEdges.empty())
        {
            if (shouldForceLocalBFS(T_init, T_end, k, S.size()))
            {
                if (shouldUseInternalBFS(T_init, k))
                {
                    int bfsDist = profiledMinRotationsBFS(T_init, T_end, k + 1);
                    if (bfsDist >= 0 && bfsDist <= k)
                        return memoReturn(true);
                }
            }

            struct FreeChoice
            {
                int parent;
                int child;
                std::pair<int,int> inserted;
            };
            std::vector<FreeChoice> choices;
            choices.reserve(freeEdges.size());
            for (const auto &edge : freeEdges)
            {
                auto inserted = insertedDiagonalAfterRotate(T_init, edge.first, edge.second);
                if (!inserted)
                    continue;
                choices.push_back(FreeChoice{edge.first, edge.second, *inserted});
            }

            std::sort(choices.begin(), choices.end(), [&](const FreeChoice &a, const FreeChoice &b) {
                const int m = static_cast<int>(T_init.original_inorder.size());
                int wa = a.inserted.second - a.inserted.first;
                int wb = b.inserted.second - b.inserted.first;
                int ba = std::min(wa, std::max(0, m - wa));
                int bb = std::min(wb, std::max(0, m - wb));
                // Prefer balanced decompositions first (shrinks both subproblems).
                if (ba != bb) return ba > bb;
                // Tie-breaker: wider inserted diagonals.
                if (wa != wb) return wa > wb;
                if (a.inserted != b.inserted) return a.inserted < b.inserted;
                if (a.parent != b.parent) return a.parent < b.parent;
                return a.child < b.child;
            });

            // Li–Xia's TreeDistS step for "free" edges is deterministic (pick a
            // free edge, rotate it, then partition). Branching over multiple
            // free edges multiplies the search tree and is the dominant source
            // of blow-ups as k grows. We keep an escape hatch for experiments,
            // but default to a single best-ranked choice.
            size_t free_edge_branch_cap = 1;
            if (const char *env = std::getenv("FLIPDIST_FREE_EDGE_BRANCH_CAP"))
            {
                int v = std::atoi(env);
                if (v > 0)
                    free_edge_branch_cap = static_cast<size_t>(v);
            }
            if (choices.size() > free_edge_branch_cap)
                choices.resize(free_edge_branch_cap);

            for (const auto &choice : choices)
            {
                if (profile::enabled())
                    ++profile::g.total_freeEdgeAttempts;
                int parent = choice.parent;
                int child = choice.child;
                debugPrint("TreeDistS: attempting free edge (" + std::to_string(parent) + "," + std::to_string(child) + ")");

                VectorRangeTreeMap T_bar = safeCopyTree(T_init);
                if (T_bar.original_nodes.empty())
                    continue;

                try
                {
                    // Rotate the candidate free edge. After a rotation, the old child becomes
                    // the parent of the edge, and the old parent becomes the child. Li–Xia’s
                    // decomposition must cut along the *new* oriented edge (parentAfter→childAfter),
                    // where childAfter realises the inserted diagonal.
                    if (T_bar.getLeftChild(parent) == child)
                    {
                        T_bar.rotateRight(parent);
                    }
                    else if (T_bar.getRightChild(parent) == child)
                    {
                        T_bar.rotateLeft(parent);
                    }
                    else
                    {
                        debugPrint("TreeDistS: Invalid rotation for candidate free edge, skipping");
                        continue;
                    }

                    int u = -1;
                    int v = -1;
                    if (!findEdgeOrientation(T_bar, parent, child, u, v))
                    {
                        debugPrint("TreeDistS: free edge endpoints not adjacent after rotation, skipping");
                        continue;
                    }

                    // Use the current lower bound as our progress proxy (do not use any
                    // externally-supplied conflict set here; it may be stale after rotations).
                    int conflicts_before = lowerBoundEdgeDifference(T_init, T_end);

                    if (TreesEqual(T_bar, T_end))
                        return memoReturn(true);

                    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S_filtered;
                    auto u_pair = twoOtherIncidentEdgesOrBoundary(T_bar, u, v);
                    if (u_pair.size() == 2)
                        appendValidPair(S_filtered, u_pair[0], u_pair[1]);

                    auto v_pair = twoOtherIncidentEdgesOrBoundary(T_bar, v, u);
                    if (v_pair.size() == 2)
                        appendValidPair(S_filtered, v_pair[0], v_pair[1]);

                    deduplicateDiagonalPairs(S_filtered, "TreeDistS_freeEdge");

                    auto parent_range = T_bar.getRange(u);
                    auto child_range = T_bar.getRange(v);
                    debugPrint("TreeDistS: partition ranges parent_range=(" + std::to_string(parent_range.first) + "," + std::to_string(parent_range.second) +
                               ") child_range=(" + std::to_string(child_range.first) + "," + std::to_string(child_range.second) + ")");

                    try
                    {
                        VectorRangeTreeMap T_bar1;
                        VectorRangeTreeMap T_bar2;
                        VectorRangeTreeMap T_end1;
                        VectorRangeTreeMap T_end2;
                        std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S1;
                        std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S2;
                        int n1 = 0;
                        int n2 = 0;

                        const std::string partitionKey = makePartitionEdgeCacheKey(T_bar, T_end, parent_range, child_range);
                        if (g_partitionFailureMemo.count(partitionKey))
                        {
                            debugPrint("TreeDistS: cached partition mismatch, skipping free edge");
                            continue;
                        }

                        auto itPart = g_partitionSuccessMemo.find(partitionKey);
                        if (itPart != g_partitionSuccessMemo.end())
                        {
                            const auto &entry = itPart->second;
                            T_bar1 = entry.startInside;
                            T_bar2 = entry.startOutside;
                            T_end1 = entry.targetInside;
                            T_end2 = entry.targetOutside;
                            S1 = entry.S_inside;
                            S2 = entry.S_outside;
                            n1 = entry.edgesInside;
                            n2 = entry.edgesOutside;
                        }
                        else
                        {
                            {
                                // Only time the actual partition construction (tree rebuilds).
                                profile::ScopedTimer partitionTimer(&profile::g.ns_partitionAlongEdge);
                                if (profile::enabled())
                                    profile::g.calls_partitionAlongEdge += 2;
                                std::tie(T_bar1, T_bar2) = VectorRangeTreeMap::partitionAlongEdge(T_bar, parent_range, child_range);
                                std::tie(T_end1, T_end2) = VectorRangeTreeMap::partitionAlongEdge(T_end, parent_range, child_range);
                            }

                            if (T_bar1.original_nodes != T_end1.original_nodes ||
                                T_bar2.original_nodes != T_end2.original_nodes)
                            {
                                if (profile::enabled())
                                    ++profile::g.total_partitionMismatches;
                                debugPrint("TreeDistS: partition mismatch, skipping free edge");
                                g_partitionFailureMemo.insert(partitionKey);
                                maybeTrimPartitionCache();
                                continue;
                            }

                            std::tie(S1, S2) = partitionS(S_filtered, T_bar1, T_bar2);
                            n1 = countInternalEdges(T_bar1);
                            n2 = countInternalEdges(T_bar2);

                            PartitionCacheEntry entry;
                            entry.startInside = T_bar1;
                            entry.startOutside = T_bar2;
                            entry.targetInside = T_end1;
                            entry.targetOutside = T_end2;
                            entry.S_inside = S1;
                            entry.S_outside = S2;
                            entry.edgesInside = n1;
                            entry.edgesOutside = n2;
                            g_partitionSuccessMemo.emplace(partitionKey, std::move(entry));
                            maybeTrimPartitionCache();
                        }
                        debugPrint("TreeDistS: partition ok n1=" + std::to_string(n1) + " n2=" + std::to_string(n2) +
                                   " |S1|=" + std::to_string(S1.size()) + " |S2|=" + std::to_string(S2.size()));

                    int conflicts_after  = lowerBoundEdgeDifference(T_bar, T_end);
                    debugPrint("TreeDistS: free-edge conflicts " + std::to_string(conflicts_before) + " -> " + std::to_string(conflicts_after));
                    // Allow plateaus; only skip if the rotation makes things worse.
                    if (conflicts_after > conflicts_before)
                    {
                        debugPrint("TreeDistS: free edge increased conflicts, skipping");
                        continue;
                    }
                    if (conflicts_after > (k - 1))
                    {
                        debugPrint("TreeDistS: free edge lb_after exceeds remaining budget, skipping");
                        continue;
                    }

                        if (n1 == 0)
                        {
                            int k2 = k - 1;
                            if (k2 >= 0 && TreeDistS(T_bar2, T_end2, k2, S2, allow_independent_retry, nullptr))
                                return memoReturn(true);
                            continue;
                        }

                        if (n2 == 0)
                        {
                            int k1 = k - 1;
                            if (k1 >= 0 && TreeDistS(T_bar1, T_end1, k1, S1, allow_independent_retry, nullptr))
                                return memoReturn(true);
                            continue;
                        }

                        // Additivity along the inserted diagonal: we only need the *minimum*
                        // budget for one side; if the other side can't be solved with the
                        // remaining budget, no split can succeed.
                        const bool solveInsideFirst = (n1 <= n2);
                        if (solveInsideFirst)
                        {
                            int k1_min = minBudgetForTreeDistS(T_bar1, T_end1, S1, k - 1, allow_independent_retry);
                            if (k1_min < 0)
                                continue;
                            int k2 = k - 1 - k1_min;
                            if (k2 >= 0 && TreeDistS(T_bar2, T_end2, k2, S2, allow_independent_retry, nullptr))
                                return memoReturn(true);
                        }
                        else
                        {
                            int k2_min = minBudgetForTreeDistS(T_bar2, T_end2, S2, k - 1, allow_independent_retry);
                            if (k2_min < 0)
                                continue;
                            int k1 = k - 1 - k2_min;
                            if (k1 >= 0 && TreeDistS(T_bar1, T_end1, k1, S1, allow_independent_retry, nullptr))
                                return memoReturn(true);
                        }
                    }
                    catch (...)
                    {
                        debugPrint("TreeDistS: exception during free-edge partition, skipping");
                        continue;
                    }
                }
                catch (...)
                {
                    debugPrint("TreeDistS: exception rotating free edge, skipping");
                    continue;
                }
            }
        }
    }

    // No budget for free-edge search or exhausted
    if (attemptedFreeEdge)
        debugPrint("TreeDistS: Free edge candidates exhausted");

    // If we have no budget, abandon this branch unless trees are already equal (handled above).
    if (k <= 0)
    {
        debugPrint("TreeDistS: No budget left");
        return memoReturn(false);
    }

    // Opportunistic local BFS only for very small subproblems.
    if (k <= 15 && shouldUseInternalBFS(T_init, k))
    {
        int bfsDist = profiledMinRotationsBFS(T_init, T_end, k);
        if (bfsDist >= 0 && bfsDist <= k)
            return memoReturn(true);
    }

    // Build a diagonal→node map for branching/rotation bookkeeping.
    auto diagMap = buildDiagonalNodeMap(T_init);

    // Conflict-focused branching: start with existing S, then add wedges for top conflicts and partner wedges that touch them.
    auto conflictsList = collectConflictingEdges(T_init, T_end);
    if (conflictsList.empty())
        return memoReturn(true);

    std::sort(conflictsList.begin(), conflictsList.end(), [](const auto &a, const auto &b) {
        int wa = a.diag.second - a.diag.first;
        int wb = b.diag.second - b.diag.first;
        if (wa != wb) return wa > wb;
        return a.diag < b.diag;
    });

    if (DEBUG && conflict_dump_budget > 0)
    {
        --conflict_dump_budget;
        std::ostringstream oss;
        oss << "TreeDistS: conflicts(" << conflictsList.size() << ") ";
        size_t cap = std::min<size_t>(5, conflictsList.size());
        for (size_t i = 0; i < cap; ++i)
        {
            const auto &c = conflictsList[i];
            oss << "diag(" << c.diag.first << "," << c.diag.second << ") edge(" << c.parent << "," << c.child << ") ";
        }
        debugPrint(oss.str());
    }

    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> branchPairs = S;

    // Keep a concrete conflict set to hand down to TreeDistI for relative comparisons
    // and to filter S-pairs. Conflicts are the diagonals present in the current
    // start tree but absent in the target.
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> conflictDiagSet;
    conflictDiagSet.reserve(conflictsList.size() * 2);
    for (const auto &c : conflictsList)
        conflictDiagSet.insert(c.diag);

    // Only synthesize extra pairs when S is empty. Once we're inside Li–Xia's
    // recursion, S already encodes the bounded neighborhood to branch on; adding
    // additional "helper" pairs can explode the branching factor at larger k.
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> extraPairs;
    size_t partnerPairsTotal = 0;
    if (S.empty())
    {
        extraPairs.reserve(64);

        auto emitWedgesForConflict = [&](const DiagonalEdge &conflict) {
            const auto conflictEdge = std::make_pair(conflict.parent, conflict.child);
            auto emitAtNode = [&](int node) {
                auto inc = getIncidentEdges(T_init, node);
                std::vector<std::pair<int,int>> others;
                others.reserve(2);
                for (const auto &e : inc)
                {
                    if ((e.first == conflictEdge.first && e.second == conflictEdge.second) ||
                        (e.first == conflictEdge.second && e.second == conflictEdge.first))
                        continue;
                    others.push_back(e);
                }
                if (others.size() < 2)
                    return;

                auto scoreEdge = [&](const std::pair<int,int> &e) -> int {
                    auto d = T_init.getRange(e.second);
                    return d.second - d.first;
                };

                if (others.size() > 2)
                {
                    std::sort(others.begin(), others.end(), [&](const auto &a, const auto &b) {
                        int wa = scoreEdge(a);
                        int wb = scoreEdge(b);
                        if (wa != wb) return wa > wb;
                        return a < b;
                    });
                    others.resize(2);
                }

                auto a = makeDiagonalEdge(T_init, others[0].first, others[0].second);
                auto b = makeDiagonalEdge(T_init, others[1].first, others[1].second);
                appendValidPair(extraPairs, a, b);
            };
            emitAtNode(conflict.parent);
            emitAtNode(conflict.child);
        };

        const size_t conflict_cap = 3;
        for (size_t i = 0; i < std::min(conflict_cap, conflictsList.size()); ++i)
            emitWedgesForConflict(conflictsList[i]);

        // Add partner wedges that touch any current conflict.
        auto partnerPairs = buildPartnerPairs(T_init, T_end);
        partnerPairsTotal = partnerPairs.size();
        for (const auto &p : partnerPairs)
        {
            if (conflictDiagSet.count(p.first.diag) || conflictDiagSet.count(p.second.diag))
            {
                extraPairs.push_back(p);
            }
        }
    }

    if (DEBUG && partner_dump_budget > 0)
    {
        --partner_dump_budget;
        std::ostringstream oss;
        oss << "TreeDistS: partnerPairs total=" << partnerPairsTotal
            << " branchPairs(after add)=" << branchPairs.size();
        size_t cap = std::min<size_t>(5, branchPairs.size());
        oss << " sample:";
        for (size_t i = 0; i < cap; ++i)
        {
            const auto &p = branchPairs[i];
            oss << " [(" << p.first.diag.first << "," << p.first.diag.second << ")/(" << p.first.parent << "," << p.first.child << ")"
                << " ; (" << p.second.diag.first << "," << p.second.diag.second << ")/(" << p.second.parent << "," << p.second.child << ")]";
        }
        debugPrint(oss.str());
    }

    deduplicateDiagonalPairs(branchPairs, "TreeDistS_conflict_base");
    deduplicateDiagonalPairs(extraPairs, "TreeDistS_conflict_extra");

    // Keep all pairs derived from S (Li–Xia), but cap the extra synthesized pairs
    // to avoid branch explosion on larger instances.
    // Cap the synthesized (non-Li–Xia) extras to avoid blow-ups, but keep
    // enough room that we don't prune away the unique wedge needed for
    // completeness at moderate k (validated via Java BFS parity).
    const size_t total_cap = std::min<size_t>(80, std::max<size_t>(16, 4ull * static_cast<size_t>(k) + 8ull));
    size_t extra_keep = 0;
    if (branchPairs.size() < total_cap)
        extra_keep = std::min(extraPairs.size(), total_cap - branchPairs.size());
    // When S is already populated (i.e., we are inside the Li–Xia recursion),
    // additional synthesized pairs can explode branching. Keep a small number
    // as a safety net, but heavily cap them for performance.
    if (!S.empty() && extra_keep > 12)
        extra_keep = 12;

    std::unordered_set<std::pair<int,int>, PairHash, PairEq> conflictsSetCurrent = conflictDiagSet;

    // Sort extras by max span (tackle widest conflicts first).
    std::sort(extraPairs.begin(), extraPairs.end(), [](const auto &a, const auto &b) {
        int wa = std::max(a.first.diag.second - a.first.diag.first,
                          a.second.diag.second - a.second.diag.first);
        int wb = std::max(b.first.diag.second - b.first.diag.first,
                          b.second.diag.second - b.second.diag.first);
        if (wa != wb) return wa > wb;
        int sa = std::min(a.first.diag.second - a.first.diag.first,
                          a.second.diag.second - a.second.diag.first);
        int sb = std::min(b.first.diag.second - b.first.diag.first,
                          b.second.diag.second - b.second.diag.first);
        return sa > sb;
    });

    if (!S.empty() && extra_keep > 8)
        extra_keep = 8;

    if (extra_keep > 0)
    {
        if (extraPairs.size() > extra_keep)
            extraPairs.resize(extra_keep);
        branchPairs.insert(branchPairs.end(), extraPairs.begin(), extraPairs.end());
        deduplicateDiagonalPairs(branchPairs, "TreeDistS_conflict_merge");
    }

    // Drop pairs that cannot contribute any conflicting edge. This is consistent
    // with the "never flip common diagonals" property: if both options are
    // already present in the target, selecting them is unnecessary.
    {
        std::vector<std::pair<DiagonalEdge, DiagonalEdge>> filtered;
        filtered.reserve(branchPairs.size());
        for (const auto &p : branchPairs)
        {
            if (conflictDiagSet.count(p.first.diag) || conflictDiagSet.count(p.second.diag))
                filtered.push_back(p);
        }
        if (filtered.size() != branchPairs.size())
            branchPairs.swap(filtered);
    }

    // Process pairs that "cover" many current conflicts first. This is purely an
    // ordering heuristic (does not prune any branch), but it has an outsized
    // impact on runtime because TreeDistS stops as soon as it finds a witness.
    std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> coverCache;
    coverCache.reserve(branchPairs.size() * 2);
    auto coverFor = [&](const std::pair<int,int> &diag) -> int {
        if (diag.second - diag.first <= 1)
            return 0;
        auto it = coverCache.find(diag);
        if (it != coverCache.end())
            return it->second;
        int count = 0;
        for (const auto &c : conflictsList)
        {
            if (c.diag.first >= diag.first && c.diag.second <= diag.second)
                ++count;
        }
        coverCache.emplace(diag, count);
        return count;
    };
    auto edgePriority = [&](const DiagonalEdge &edge) -> long long {
        int cover = coverFor(edge.diag);
        int width = edge.diag.second - edge.diag.first;
        // Weight coverage heavily; width is a secondary tie-breaker.
        return static_cast<long long>(cover) * 10'000ll + static_cast<long long>(width);
    };
    auto pairPriority = [&](const std::pair<DiagonalEdge, DiagonalEdge> &p) {
        int c1 = coverFor(p.first.diag);
        int c2 = coverFor(p.second.diag);
        int w1 = p.first.diag.second - p.first.diag.first;
        int w2 = p.second.diag.second - p.second.diag.first;
        int cover_union = c1 + c2 - std::min(c1, c2) / 2; // discount overlap slightly
        int cover_max = std::max(c1, c2);
        int width_max = std::max(w1, w2);
        int width_min = std::min(w1, w2);
        // Prefer pairs that meaningfully cover conflicts; zero-coverage pairs
        // should be deprioritised (though most were already filtered).
        return std::make_tuple(cover_union, cover_max, width_max, width_min);
    };
    for (auto &pair : branchPairs)
    {
        if (edgePriority(pair.second) > edgePriority(pair.first))
            std::swap(pair.first, pair.second);
    }
    std::sort(branchPairs.begin(), branchPairs.end(), [&](const auto &a, const auto &b) {
        auto ta = pairPriority(a);
        auto tb = pairPriority(b);
        if (ta != tb) return ta > tb;
        long long a1 = edgePriority(a.first);
        long long b1 = edgePriority(b.first);
        if (a1 != b1) return a1 > b1;
        long long a2 = edgePriority(a.second);
        long long b2 = edgePriority(b.second);
        if (a2 != b2) return a2 > b2;
        // tie-breaker: deterministic key
        return makeDiagonalPairKey(a.first, a.second) < makeDiagonalPairKey(b.first, b.second);
    });

    if (branchPairs.empty())
        return memoReturn(false);

    bool result = false;
    const int maxNode = std::max(T_init.max_node_value, T_end.max_node_value);
    if (maxNode > 0 && maxNode < 64)
    {
        // Branch memoization must distinguish different chosen edge-sets. Build a stable
        // ID map for the oriented parent->child edges that can be selected from S.
        std::vector<std::pair<int,int>> selectableEdges;
        selectableEdges.reserve(branchPairs.size() * 2);
        for (const auto &pair : branchPairs)
        {
            const DiagonalEdge opts[2] = {pair.first, pair.second};
            for (const auto &opt : opts)
            {
                auto oriented = orientDiagonalEdge(T_init, opt, diagMap);
                if (!oriented)
                {
                    if (hasParentChildEdge(T_init, opt.parent, opt.child))
                        oriented = std::make_pair(opt.parent, opt.child);
                }
                if (!oriented)
                    continue;
                if (!hasParentChildEdge(T_init, oriented->first, oriented->second))
                    continue;
                selectableEdges.push_back(*oriented);
            }
        }
        std::sort(selectableEdges.begin(), selectableEdges.end());
        selectableEdges.erase(std::unique(selectableEdges.begin(), selectableEdges.end()), selectableEdges.end());

        std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> edgeIdMap;
        edgeIdMap.reserve(selectableEdges.size());
        for (size_t i = 0; i < selectableEdges.size(); ++i)
        {
            edgeIdMap.emplace(selectableEdges[i], static_cast<int>(i));
        }

        std::unordered_set<BranchMemoKeyMask, BranchMemoKeyMaskHash> branchMemo;
        branchMemo.reserve(branchPairs.size() * 8);
        std::vector<std::pair<int,int>> chosen;
        chosen.reserve(branchPairs.size());
        result = branchOnSPairsMask(T_init, T_end, k, branchPairs, 0, chosen, 0ull, 0ull, 0ull,
                                    branchMemo, diagMap, edgeIdMap, allow_independent_retry, &conflictsSetCurrent);
    }
    else
    {
        std::unordered_set<std::string> branchMemo;
        std::vector<std::pair<int,int>> chosen;
        chosen.reserve(branchPairs.size());
        std::unordered_set<int> used_nodes;
        result = branchOnSPairsSet(T_init, T_end, k, branchPairs, 0, chosen, used_nodes, branchMemo, diagMap, allow_independent_retry, &conflictsSetCurrent);
    }
    if (!result && k > 0)
    {
        // Last-resort greedy step: rotate the realizing edge of the widest conflict and recurse with k-1.
        auto conflicts = collectConflictingEdges(T_init, T_end);
        if (!conflicts.empty())
        {
            auto widest = *std::max_element(conflicts.begin(), conflicts.end(), [](const auto &a, const auto &b) {
                int wa = a.diag.second - a.diag.first;
                int wb = b.diag.second - b.diag.first;
                return wa < wb;
            });
            int parent = widest.parent, child = widest.child;
            VectorRangeTreeMap T_greedy = safeCopyTree(T_init);
            if (T_greedy.original_nodes.size() > 0)
            {
                if (!hasParentChildEdge(T_greedy, parent, child))
                {
                    int p, c;
                    if (findEdgeOrientation(T_greedy, parent, child, p, c))
                    {
                        parent = p;
                        child = c;
                    }
                }
                if (hasParentChildEdge(T_greedy, parent, child))
                {
                    if (T_greedy.getLeftChild(parent) == child)
                        T_greedy.rotateRight(parent);
                    else if (T_greedy.getRightChild(parent) == child)
                        T_greedy.rotateLeft(parent);
                    int conf_before = lowerBoundEdgeDifference(T_init, T_end);
                    int conf_after = lowerBoundEdgeDifference(T_greedy, T_end);
                    if (conf_after < conf_before)
                    {
                        if (TreeDistS(T_greedy, T_end, k - 1, {}, allow_independent_retry, &conflictsSetCurrent))
                            result = true;
                    }
                }
            }
        }
    }
    return memoReturn(result);
}


} // namespace flipdist
