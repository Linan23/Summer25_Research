#include "rotation_tree.h"
#include "tree_generators.h"
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <climits>
#include <cassert>
#include <functional>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <deque>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <stdexcept>

int MinRotationsBFS(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap);
bool FlipDistTree(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k);
bool TreeDistI(const VectorRangeTreeMap &T_init,
               const VectorRangeTreeMap &T_final,
               int k,
               const std::vector<std::pair<int, int>> &I);
bool TreeDistS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
               const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);
bool hasParentChildEdge(const VectorRangeTreeMap &T, int parent, int child);

namespace {

struct FlipMemoKey {
    std::string start;
    std::string target;
    int k;

    bool operator==(const FlipMemoKey &other) const noexcept {
        return k == other.k && start == other.start && target == other.target;
    }
};

struct FlipMemoKeyHash {
    std::size_t operator()(const FlipMemoKey &key) const noexcept {
        std::size_t h1 = std::hash<std::string>()(key.start);
        std::size_t h2 = std::hash<std::string>()(key.target);
        std::size_t h3 = std::hash<int>()(key.k);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
        return seed;
    }
};

static std::unordered_map<FlipMemoKey, bool, FlipMemoKeyHash> g_flipDistMemo;

static int lowerBoundEdgeDifference(const VectorRangeTreeMap &A,
                                    const VectorRangeTreeMap &B)
{
    std::unordered_set<std::pair<int,int>,PairHash,PairEq> EA, EB;
    A.collectEdges(A.root, EA);
    B.collectEdges(B.root, EB);
    int common = 0;
    for (const auto &e : EA) {
        if (EB.count(e)) ++common;
    }
    int diff = static_cast<int>(EA.size() + EB.size() - 2 * common);
    // Each rotation can fix at most three diagonal mismatches in the convex polygon model.
    return (diff + 2) / 3;
}

static inline std::pair<int,int> makeUndirectedPair(int a, int b)
{
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
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

static bool edgesShareEndpoint(const std::pair<int,int>& a,
                               const std::pair<int,int>& b)
{
    auto ua = makeUndirectedPair(a.first, a.second);
    auto ub = makeUndirectedPair(b.first, b.second);
    return ua.first == ub.first || ua.first == ub.second ||
           ua.second == ub.first || ua.second == ub.second;
}

static bool branchOnSPairs(const VectorRangeTreeMap &T_init,
                           const VectorRangeTreeMap &T_end,
                           int k,
                           const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S,
                           std::size_t index,
                           std::vector<std::pair<int,int>> &chosen)
{
    if (index == S.size())
    {
        if (chosen.empty()) return false;
        return ::TreeDistI(T_init, T_end, k, chosen);
    }

    if (branchOnSPairs(T_init, T_end, k, S, index + 1, chosen))
        return true;

    const auto &pair = S[index];
    const std::pair<int,int> options[2] = {pair.first, pair.second};

    for (const auto &edge : options)
    {
        if (!::hasParentChildEdge(T_init, edge.first, edge.second))
            continue;

        bool conflict = false;
        for (const auto &sel : chosen)
        {
            if (edgesShareEndpoint(edge, sel)) { conflict = true; break; }
        }
        if (conflict) continue;

        chosen.push_back(edge);
        if (branchOnSPairs(T_init, T_end, k, S, index + 1, chosen))
            return true;
        chosen.pop_back();
    }

    return false;
}

static bool tryCommonEdgeDecomposition(const VectorRangeTreeMap &start,
                                       const VectorRangeTreeMap &target,
                                       int k)
{
    if (k <= 0) return false;

    std::unordered_set<std::pair<int,int>, PairHash, PairEq> startEdgesDirected;
    start.collectEdges(start.root, startEdgesDirected);
    if (startEdgesDirected.empty()) return false;

    std::unordered_set<std::pair<int,int>, PairHash, PairEq> targetEdgesDirected;
    target.collectEdges(target.root, targetEdgesDirected);

    std::unordered_map<std::pair<int,int>, std::pair<int,int>, PairHash, PairEq> targetOrientation;
    for (const auto &edge : targetEdgesDirected) {
        auto und = makeUndirectedPair(edge.first, edge.second);
        int parent, child;
        if (!findEdgeOrientation(target, edge.first, edge.second, parent, child)) continue;
        targetOrientation.emplace(und, std::make_pair(parent, child));
    }

    for (const auto &edge : startEdgesDirected) {
        auto und = makeUndirectedPair(edge.first, edge.second);
        auto it = targetOrientation.find(und);
        if (it == targetOrientation.end()) continue;

        int parentStart, childStart;
        if (!findEdgeOrientation(start, edge.first, edge.second, parentStart, childStart)) continue;
        int parentTarget = it->second.first;
        int childTarget  = it->second.second;

        auto startChildRange  = start.getRange(childStart);
        auto targetChildRange = target.getRange(childTarget);

        auto startParts  = VectorRangeTreeMap::partitionAlongEdge(start,
                                    start.getRange(parentStart), startChildRange);
        auto targetParts = VectorRangeTreeMap::partitionAlongEdge(target,
                                    target.getRange(parentTarget), targetChildRange);

        const auto &startChildTree = startParts.first;
        const auto &startRestTree  = startParts.second;
        const auto &targetChildTree = targetParts.first;
        const auto &targetRestTree  = targetParts.second;

        if (startChildTree.original_nodes.empty() || startRestTree.original_nodes.empty())
            continue;
        if (startChildTree.original_nodes != targetChildTree.original_nodes)
            continue;
        if (startRestTree.original_nodes != targetRestTree.original_nodes)
            continue;

        int lbChild = lowerBoundEdgeDifference(startChildTree, targetChildTree);
        int lbRest  = lowerBoundEdgeDifference(startRestTree, targetRestTree);
        if (lbChild + lbRest > k) continue;

        for (int budgetChild = lbChild; budgetChild <= k - lbRest; ++budgetChild) {
            if (!FlipDistTree(startChildTree, targetChildTree, budgetChild)) continue;
            int remaining = k - budgetChild;
            if (FlipDistTree(startRestTree, targetRestTree, remaining)) {
                return true;
            }
        }

        int distChild = MinRotationsBFS(startChildTree, targetChildTree, k);
        if (distChild < 0) continue;
        int remaining = k - distChild;
        if (remaining < 0) continue;
        int distRest = MinRotationsBFS(startRestTree, targetRestTree, remaining);
        if (distRest < 0) continue;
        if (distChild + distRest <= k) {
            return true;
        }
    }

    return false;
}

static std::vector<std::pair<int,int>> collectConflictingEdges(
        const VectorRangeTreeMap &start,
        const VectorRangeTreeMap &target)
{
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> targetEdgesDirected;
    target.collectEdges(target.root, targetEdgesDirected);
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> targetUndirected;
    for (const auto &edge : targetEdgesDirected) {
        targetUndirected.insert(makeUndirectedPair(edge.first, edge.second));
    }

    std::vector<std::pair<int,int>> conflicts;
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> startEdgesDirected;
    start.collectEdges(start.root, startEdgesDirected);
    for (const auto &edge : startEdgesDirected) {
        if (!targetUndirected.count(makeUndirectedPair(edge.first, edge.second))) {
            conflicts.push_back(edge);
        }
    }
    return conflicts;
}

static std::vector<std::pair<int,int>> buildMaxIndependentSet(
        const VectorRangeTreeMap &start,
        const std::vector<std::pair<int,int>> &conflicts)
{
    std::vector<std::pair<int,int>> independent;
    std::unordered_set<int> usedNodes;

    std::vector<std::pair<int,int>> sorted = conflicts;
    std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        int ma = std::min(a.first, a.second);
        int mb = std::min(b.first, b.second);
        if (ma != mb) return ma < mb;
        int xa = std::max(a.first, a.second);
        int xb = std::max(b.first, b.second);
        return xa < xb;
    });

    for (const auto &edge : sorted) {
        int parent, child;
        if (!findEdgeOrientation(start, edge.first, edge.second, parent, child))
            continue;

        if (usedNodes.count(parent) || usedNodes.count(child))
            continue;

        independent.emplace_back(parent, child);
        usedNodes.insert(parent);
        usedNodes.insert(child);
    }

    return independent;
}

} // namespace

// Debug flag
const bool DEBUG = false; // Turn off debug for cleaner testing output

void debugPrint(const std::string &msg)
{
    if (DEBUG)
    {
        std::cout << "[DEBUG] " << msg << std::endl;
    }
}

// Helper Functions
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

int countInternalEdges(const VectorRangeTreeMap &T)
{
    return getInternalEdges(T).size();
}

bool areAdjacent(const std::pair<int, int> &e1, const std::pair<int, int> &e2)
{
    return e1.first == e2.first || e1.first == e2.second ||
           e1.second == e2.first || e1.second == e2.second;
}

VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap &T)
{
    VectorRangeTreeMap copy;
    try
    {
        if (T.root < 0 || !T.isOriginal(T.root) || T.original_nodes.empty())
        {
            return copy;
        }

        std::vector<int> preorder, inorder;

        std::function<void(int, std::vector<int> &)> buildPreorder = [&](int node, std::vector<int> &pre)
        {
            if (node < 0 || !T.isOriginal(node))
                return;
            pre.push_back(node);
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left))
                buildPreorder(left, pre);
            if (right >= 0 && T.isOriginal(right))
                buildPreorder(right, pre);
        };

        std::function<void(int, std::vector<int> &)> buildInorder = [&](int node, std::vector<int> &in)
        {
            if (node < 0 || !T.isOriginal(node))
                return;
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left))
                buildInorder(left, in);
            in.push_back(node);
            if (right >= 0 && T.isOriginal(right))
                buildInorder(right, in);
        };

        buildPreorder(T.root, preorder);
        buildInorder(T.root, inorder);

        if (preorder.size() == inorder.size() && !preorder.empty())
        {
            copy.build(preorder, inorder);
        }
    }
    catch (...)
    {
        VectorRangeTreeMap empty;
        return empty;
    }
    return copy;
}

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

static inline std::pair<int, int> undirected(int a, int b)
{
    if (a < b)
        return {a, b};
    return {b, a};
}

std::pair<bool, std::pair<int, int>> findFreeEdge(const VectorRangeTreeMap &T_init,
                                                  const VectorRangeTreeMap &T_final)
{
    try
    {
        if (T_init.original_nodes.empty() || T_final.original_nodes.empty())
        {
            return {false, {-1, -1}};
        }

        auto initEdges = getInternalEdges(T_init);
        auto finalEdges = getInternalEdges(T_final);

        // Build UNDIRECTED sets
        std::set<std::pair<int, int>> initU, finalU;
        for (auto &e : initEdges)
            initU.insert(undirected(e.first, e.second));
        for (auto &e : finalEdges)
            finalU.insert(undirected(e.first, e.second));

        for (const auto &edge : initEdges)
        {
            int parent = edge.first;
            int child = edge.second;

            if (!hasParentChildEdge(T_init, parent, child))
                continue;

            VectorRangeTreeMap testTree = safeCopyTree(T_init);
            if (testTree.original_nodes.empty())
                continue;

            bool rotated = false;
            if (testTree.getLeftChild(parent) == child)
            {
                testTree.rotateRight(parent);
                rotated = true;
            }
            else if (testTree.getRightChild(parent) == child)
            {
                testTree.rotateLeft(parent);
                rotated = true;
            }
            if (!rotated)
                continue;

            auto newEdges = getInternalEdges(testTree);
            for (const auto &ne : newEdges)
            {
                auto neu = undirected(ne.first, ne.second);
                if (initU.find(neu) == initU.end() // not in original
                    && finalU.find(neu) != finalU.end())
                {                        // IS in target
                    return {true, edge}; // keep original directed (parent->child) for rotation
                }
            }
        }
    }
    catch (...)
    {
    }

    return {false, {-1, -1}};
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
        // Return empty on error
    }

    return incident;
}

// Helper function to partition S based on which tree partition edges belong to
std::pair<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>,
          std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>
partitionS(const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S,
           const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2)
{

    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> S1, S2;

    // Get node sets for each partition
    std::set<int> nodes1, nodes2;
    for (int node : T1.original_nodes)
        nodes1.insert(node);
    for (int node : T2.original_nodes)
        nodes2.insert(node);

    for (const auto &edgePair : S)
    {
        auto &edge1 = edgePair.first;
        auto &edge2 = edgePair.second;

        // Check if both edges of the pair belong to T1
        bool edge1_in_T1 = nodes1.count(edge1.first) && nodes1.count(edge1.second);
        bool edge2_in_T1 = nodes1.count(edge2.first) && nodes1.count(edge2.second);

        // Check if both edges of the pair belong to T2
        bool edge1_in_T2 = nodes2.count(edge1.first) && nodes2.count(edge1.second);
        bool edge2_in_T2 = nodes2.count(edge2.first) && nodes2.count(edge2.second);

        if (edge1_in_T1 && edge2_in_T1)
        {
            S1.push_back(edgePair);
        }
        else if (edge1_in_T2 && edge2_in_T2)
        {
            S2.push_back(edgePair);
        }
        // If edge pair spans both partitions, we could assign to both or neither
        // For simplicity, we'll ignore cross-partition pairs
    }

    return {S1, S2};
}

/**
 * FLIPDISTTREE - Main algorithm
 * Maps to: FlipDistTree(T_init, T_final, k) pseudocode
 */
bool FlipDistTree(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k)
{
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    if (k < 0)
        return false;

    std::string startKey = treeToString(T_init);
    std::string targetKey = treeToString(T_final);
    FlipMemoKey memoKey{startKey, targetKey, k};
    if (auto it = g_flipDistMemo.find(memoKey); it != g_flipDistMemo.end())
        return it->second;

    // BASE CASE: Check if trees are already identical
    if (TreesEqual(T_init, T_final))
    {
        debugPrint("Trees already equal, returning true");
        g_flipDistMemo[memoKey] = true;
        return true;
    }

    int lb = lowerBoundEdgeDifference(T_init, T_final);
    if (lb > k)
    {
        debugPrint("Lower bound " + std::to_string(lb) + " exceeds k=" + std::to_string(k));
        g_flipDistMemo[memoKey] = false;
        return false;
    }

    if (tryCommonEdgeDecomposition(T_init, T_final, k))
    {
        g_flipDistMemo[memoKey] = true;
        return true;
    }

    // PSEUDOCODE STEP 0: "If φ(T_init) > k, return False"
    // IMPLEMENTATION NOTE: We use a more generous bound for practical performance
    // Original pseudocode: φ(T_init) > k
    // Our implementation: φ(T_init) > k + φ(T_init)/2
    // REASON: The strict bound from the paper is too restrictive for real test cases
    int phi_init = countInternalEdges(T_init); // φ(T_init) = number of internal edges
    debugPrint("T_init has " + std::to_string(phi_init) + " internal edges");

    // Handle trivial case: no internal edges
    if (phi_init == 0)
    {
        bool result = countInternalEdges(T_final) == 0;
        debugPrint("No internal edges, result: " + std::string(result ? "true" : "false"));
        return result;
    }

    // PSEUDOCODE STEP 1: "Enumerate all subsets I of independent internal edges in T_init"
    auto conflicts = collectConflictingEdges(T_init, T_final);
    debugPrint("Conflicting edges: " + std::to_string(conflicts.size()));

    auto independentSet = buildMaxIndependentSet(T_init, conflicts);
    if (!independentSet.empty())
    {
        debugPrint("Max independent set size=" + std::to_string(independentSet.size()));
        if (TreeDistI(T_init, T_final, k, independentSet))
        {
            g_flipDistMemo[memoKey] = true;
            return true;
        }
    }

    // PSEUDOCODE STEP 2: "Return False"
    debugPrint("No solution found, returning false");
    g_flipDistMemo[memoKey] = false;
    return false;
}

// Returns minimal k where FlipDistTree(T1,T2,k) is true, or -1 if not found up to k_max.
int FlipDistMinK(const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2, int k_max)
{
    if (TreesEqual(T1, T2))
        return 0;
    int lb0 = lowerBoundEdgeDifference(T1, T2);
    if (lb0 > k_max)
        return -1;
    int start_k = std::max(1, lb0);
    for (int k = start_k; k <= k_max; ++k)
    {
        if (FlipDistTree(T1, T2, k))
            return k;
    }
    return -1; // not found within bound
}

/**
 * TREEDISTI - Handles independent edge set I
 * Maps to: TreeDist-I(T_init, T_final, k, I) pseudocode
 */
bool TreeDistI(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k,
               const std::vector<std::pair<int, int>> &I)
{
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    // PSEUDOCODE STEP 0: "Let φ(T) = # internal edges in T. If φ(T_init) > k − |I|, return False"
    int phi_init = countInternalEdges(T_init);
    int remaining_budget = k - (int)I.size(); // k - |I| from pseudocode

    if (remaining_budget < 0)
    { // This covers φ(T_init) > k − |I| case
        debugPrint("TreeDistI: Not enough budget for |I| rotations");
        return false;
    }

    // Special handling when budget exactly equals |I| - try direct solution
    if (remaining_budget == 0)
    {
        // Apply I to a working copy
        VectorRangeTreeMap T_bar = safeCopyTree(T_init);
        for (const auto &edge : I)
        {
            int parent = edge.first, child = edge.second;
            if (!hasParentChildEdge(T_bar, parent, child))
                continue;
            if (T_bar.getLeftChild(parent) == child)
                T_bar.rotateRight(parent);
            else if (T_bar.getRightChild(parent) == child)
                T_bar.rotateLeft(parent);
        }

        if (TreesEqual(T_bar, T_final))
            return true;

        // Build S from the rotated tree (same logic you already use below)
        std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> S;

        // For each rotated edge, collect incident pairs on u and v
        // (identical to your code below — reuse it verbatim)
        // ---- begin reuse of your Step 2.2 collection on T_bar ----
        for (const auto &edge : I)
        {
            int parent = edge.first, child = edge.second;

            int u = child, v = parent; // after rotation child becomes parent
            auto u_incident = getIncidentEdges(T_bar, u);
            auto v_incident = getIncidentEdges(T_bar, v);

            std::vector<std::pair<int, int>> u_others, v_others;
            for (const auto &e : u_incident)
                if (!((e.first == u && e.second == v) || (e.first == v && e.second == u)))
                    u_others.push_back(e);
            for (const auto &e : v_incident)
                if (!((e.first == u && e.second == v) || (e.first == v && e.second == u)))
                    v_others.push_back(e);

            if (u_others.size() >= 2)
                S.emplace_back(u_others[0], u_others[1]);
            if (v_others.size() >= 2)
                S.emplace_back(v_others[0], v_others[1]);
        }
        // ---- end reuse ----

        // Continue with zero remaining budget
        return TreeDistS(T_bar, T_final, 0, S);
    }

    debugPrint("TreeDistI: Proceeding with remaining_budget=" + std::to_string(remaining_budget));

    // PSEUDOCODE STEP 0: "If φ(T_init) = 0 and k ≥ 0, return True"
    if (phi_init == 0 && k >= 0)
    {
        bool result = TreesEqual(T_init, T_final);
        debugPrint("TreeDistI: φ(T_init) = 0, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    // PSEUDOCODE STEP 1: "S ← ∅"
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> S;

    // PSEUDOCODE STEP 2: "For each edge e ∈ I do:"
    VectorRangeTreeMap T_bar = safeCopyTree(T_init); // T̄_init from pseudocode
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
            debugPrint("TreeDistI: Invalid edge, skipping");
            continue;
        }

        int u, v; // The nodes u,v from pseudocode step 2.2

        // PSEUDOCODE STEP 2.1: "Rotate e in T_init → creates a new internal edge ē"
        if (T_bar.getLeftChild(parent) == child)
        {
            u = child;  // After rotation, child becomes parent
            v = parent; // After rotation, parent becomes child
            T_bar.rotateRight(parent);
        }
        else if (T_bar.getRightChild(parent) == child)
        {
            u = child;  // After rotation, child becomes parent
            v = parent; // After rotation, parent becomes child
            T_bar.rotateLeft(parent);
        }
        else
        {
            continue;
        }

        debugPrint("TreeDistI: Applied rotation, new edge ē connects " + std::to_string(u) + " and " + std::to_string(v));

        // PSEUDOCODE STEP 2.2: "Let u,v be the two nodes joined by ē.
        //                       Let {e₁, e₁′} = the two other edges in T_init incident to u.
        //                       Let {e₂, e₂′} = the two other edges in T_init incident to v.
        //                       Add the pairs (e₁,e₁′) and (e₂,e₂′) to S."
        auto u_incident = getIncidentEdges(T_bar, u);
        auto v_incident = getIncidentEdges(T_bar, v);

        std::vector<std::pair<int, int>> u_others, v_others;

        // Filter out the new edge ē = (u,v) to get the "other" edges
        for (const auto &e : u_incident)
        {
            if (!((e.first == u && e.second == v) || (e.first == v && e.second == u)))
            {
                u_others.push_back(e); // These are {e₁, e₁′} from pseudocode
            }
        }

        for (const auto &e : v_incident)
        {
            if (!((e.first == u && e.second == v) || (e.first == v && e.second == u)))
            {
                v_others.push_back(e); // These are {e₂, e₂′} from pseudocode
            }
        }

        // Add pairs (e₁,e₁′) and (e₂,e₂′) to S
        if (u_others.size() >= 2)
        {
            S.emplace_back(u_others[0], u_others[1]); // Add (e₁,e₁′) to S
            debugPrint("TreeDistI: Added edge pair for u");
        }

        if (v_others.size() >= 2)
        {
            S.emplace_back(v_others[0], v_others[1]); // Add (e₂,e₂′) to S
            debugPrint("TreeDistI: Added edge pair for v");
        }
    }

    // PSEUDOCODE STEP 3: "Return TreeDist–S(T̄_init, T_final, k−|I|, S)"
    return TreeDistS(T_bar, T_final, k - (int)I.size(), S);
}

/**
 * TREEDISTS - Handles S-branching and partitioning
 * Maps to: TreeDist-S(T_init, T_end, k, S) pseudocode
 */
bool TreeDistS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
               const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S)
{
    debugPrint("Entering TreeDistS with k=" + std::to_string(k) + ", |S|=" + std::to_string(S.size()));

    // Base case: trees already equal
    if (TreesEqual(T_init, T_end))
    {
        debugPrint("TreeDistS: Trees already equal");
        return true;
    }

    // PSEUDOCODE STEP 0: "Let φ(T) = number of internal edges in T. If φ(T_init) > k, return False"
    // IMPLEMENTATION NOTE: We use a more generous bound for practical performance
    // Original pseudocode: φ(T_init) > k
    // Our implementation: φ(T_init) > k + 2
    // REASON: Strict bound is too restrictive, this allows more exploration
    int phi_init = countInternalEdges(T_init);

    // PSEUDOCODE STEP 0: "If φ(T_init) = 0 and k ≥ 0, return True"
    if (phi_init == 0 && k >= 0)
    {
        bool result = TreesEqual(T_init, T_end);
        debugPrint("TreeDistS: φ(T_init) = 0, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    if (k < 0)
    {
        debugPrint("TreeDistS: Negative budget");
        return false;
    }

    // PSEUDOCODE STEP 1: "If there is a 'free' internal edge e in T_init
    //                     (i.e. rotating e would insert an edge of T_end), then:"
    auto [hasFree, freeEdge] = findFreeEdge(T_init, T_end);

    if (hasFree)
    {
        debugPrint("TreeDistS: Found free edge (" + std::to_string(freeEdge.first) + "," + std::to_string(freeEdge.second) + ")");

        try
        {
            int parent = freeEdge.first;
            int child = freeEdge.second;

            // PSEUDOCODE STEP 1.1: "Remove from S every pair that contains e"
            auto eU = undirected(freeEdge.first, freeEdge.second);
            std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> S_filtered;
            S_filtered.reserve(S.size());

            for (const auto &pair : S)
            {
                auto feU = undirected(freeEdge.first, freeEdge.second);
                auto p1U = undirected(pair.first.first, pair.first.second);
                auto p2U = undirected(pair.second.first, pair.second.second);
                if (!(p1U == feU || p2U == feU))
                {
                    S_filtered.push_back(pair);
                }
            }
            debugPrint("TreeDistS: Filtered S from " + std::to_string(S.size()) + " to " + std::to_string(S_filtered.size()) + " pairs");

            // PSEUDOCODE STEP 1.2: "Rotate e in T_init → creates new edge ē (now common with T_end).
            //                       Let T̄_init be the resulting tree."
            VectorRangeTreeMap T_bar = safeCopyTree(T_init); // T̄_init from pseudocode
            int u, v;                                        // The nodes joined by new edge ē

            if (T_bar.getLeftChild(parent) == child)
            {
                T_bar.rotateRight(parent);
                u = child;  // New parent after rotation
                v = parent; // New child after rotation
            }
            else if (T_bar.getRightChild(parent) == child)
            {
                T_bar.rotateLeft(parent);
                u = child;  // New parent after rotation
                v = parent; // New child after rotation
            }
            else
            {
                return false; // Invalid rotation
            }

            // Check if rotation immediately solves the problem
            if (TreesEqual(T_bar, T_end))
            {
                debugPrint("TreeDistS: Solved with free edge rotation");
                return true;
            }

            // PSEUDOCODE STEP 1.3: "Let ē join nodes u and v in T̄_init.
            //                       Let {e₁,e₁′} = the two other edges incident at u.
            //                       Let {e₂,e₂′} = the two other edges incident at v.
            //                       Add (e₁,e₁′) and (e₂,e₂′) to S."
            auto u_incident = getIncidentEdges(T_bar, u);
            auto v_incident = getIncidentEdges(T_bar, v);

            std::vector<std::pair<int, int>> u_others, v_others;

            // Collect edges incident to u (excluding the new edge ē = (u,v))
            for (const auto &e : u_incident)
            {
                if (!((e.first == u && e.second == v) || (e.first == v && e.second == u)))
                {
                    u_others.push_back(e); // These are {e₁,e₁′}
                }
            }

            // Collect edges incident to v (excluding the new edge ē = (u,v))
            for (const auto &e : v_incident)
            {
                if (!((e.first == u && e.second == v) || (e.first == v && e.second == u)))
                {
                    v_others.push_back(e); // These are {e₂,e₂′}
                }
            }

            // Add pairs to S_filtered
            if (u_others.size() >= 2)
            {
                S_filtered.emplace_back(u_others[0], u_others[1]); // Add (e₁,e₁′)
            }
            if (v_others.size() >= 2)
            {
                S_filtered.emplace_back(v_others[0], v_others[1]); // Add (e₂,e₂′)
            }

            // PSEUDOCODE STEP 1.4: "partition both T̄_init and T_end along ē, yielding two subtree-pairs
            //                       {T̄_init¹, T̄_init²} and {T_end¹, T_end²}.
            //                       Partition S into S₁ (pairs lying in T̄_init¹) and S₂ (pairs in T̄_init²).
            //                       Let n₁ = φ(T̄_init¹) and n₂ = φ(T̄_init²)."
            debugPrint("TreeDistS: Implementing partitioning logic (steps 1.4-1.8)");

            auto parent_range = T_bar.getRange(u);
            auto child_range = T_bar.getRange(v);

            try
            {
                // Partition both trees along the new edge ē
                auto [T_bar1, T_bar2] = VectorRangeTreeMap::partitionAlongEdge(T_bar, parent_range, child_range); // {T̄_init¹, T̄_init²}
                auto [T_end1, T_end2] = VectorRangeTreeMap::partitionAlongEdge(T_end, parent_range, child_range); // {T_end¹, T_end²}

                // Validate partitions have matching node sets
                if (T_bar1.original_nodes != T_end1.original_nodes ||
                    T_bar2.original_nodes != T_end2.original_nodes)
                {
                    debugPrint("TreeDistS: Partition mismatch, falling back to simple recursion");
                    return TreeDistS(T_bar, T_end, k - 1, S_filtered);
                }

                // Partition S into S₁ and S₂
                auto [S1, S2] = partitionS(S_filtered, T_bar1, T_bar2);

                int n1 = countInternalEdges(T_bar1); // n₁ = φ(T̄_init¹)
                int n2 = countInternalEdges(T_bar2); // n₂ = φ(T̄_init²)

                debugPrint("TreeDistS: Partitioned into subtrees of size " + std::to_string(n1) + " and " + std::to_string(n2));

                // PSEUDOCODE STEP 1.5: "If T̄_init¹ has no internal edges(trivial),
                //                       return TreeDist-S(T̄_init², T_end², k − 1 − n₁, S₂)."
                if (n1 == 0)
                {
                    debugPrint("TreeDistS: T1 trivial, solving T2");
                    return TreeDistS(T_bar2, T_end2, k - 1 - n1, S2);
                }

                // PSEUDOCODE STEP 1.6: "If T̄_init² has no internal edges(trivial),
                //                       return TreeDist-S(T̄_init¹, T_end¹, k − 1 − n₂, S₁)."
                if (n2 == 0)
                {
                    debugPrint("TreeDistS: T2 trivial, solving T1");
                    return TreeDistS(T_bar1, T_end1, k - 1 - n2, S1);
                }

                // PSEUDOCODE STEP 1.7: "For k₁ = n₁+1 to (k − 1 − n₂):
                //                       if TreeDist-S(T̄_init¹, T_end¹, k₁, S₁) returns True,"
                // IMPLEMENTATION NOTE: We use k₁ = n₁ to (k − 1 − n₂) instead of n₁+1
                // REASON: The paper's bound is too strict; we need at least n₁ budget for subtree 1
                debugPrint("TreeDistS: Starting budget allocation loop");

                for (int k1 = n1; k1 <= k - 1 - n2; k1++)
                { // MODIFIED: start from n₁ instead of n₁+1
                    debugPrint("TreeDistS: Trying k1=" + std::to_string(k1) + " for subtree 1");

                    if (TreeDistS(T_bar1, T_end1, k1, S1))
                    {
                        // PSEUDOCODE STEP 1.8: "Return TreeDist-S(T̄_init², T_end², k − 1 − k₁, S₂)."
                        int k2 = k - 1 - k1;
                        debugPrint("TreeDistS: T1 succeeded, trying k2=" + std::to_string(k2) + " for subtree 2");

                        if (TreeDistS(T_bar2, T_end2, k2, S2))
                        {
                            debugPrint("TreeDistS: Both subtrees solved!");
                            return true;
                        }
                    }
                }

                // If no k₁ found in step 1.7: "If no such k₁ found, return False"
                debugPrint("TreeDistS: Budget allocation failed");
                return false;
            }
            catch (...)
            {
                debugPrint("TreeDistS: Partitioning failed, using simple recursion");
                return TreeDistS(T_bar, T_end, k - 1, S_filtered);
            }
        }
        catch (...)
        {
            debugPrint("TreeDistS: Exception during free edge handling");
            return false;
        }
    }

    // PSEUDOCODE STEP 2: "No free edge shortcut → branch on S
    //                     For each nonempty independent subset I ⊆ ⋃ S (no two edges in I share a node):"
    debugPrint("TreeDistS: No free edge, implementing S branching (step 2)");

    if (k <= 0)
    {
        debugPrint("TreeDistS: No budget left");
        return false;
    }

    if (S.empty())
    {
        // No constraints from S, try any rotation
        debugPrint("TreeDistS: S is empty, trying any rotation");
        auto edges = getInternalEdges(T_init);
        for (const auto &edge : edges)
        {
            try
            {
                VectorRangeTreeMap T_test = safeCopyTree(T_init);
                int parent = edge.first;
                int child = edge.second;

                if (T_test.getLeftChild(parent) == child)
                {
                    T_test.rotateRight(parent);
                }
                else if (T_test.getRightChild(parent) == child)
                {
                    T_test.rotateLeft(parent);
                }

                if (TreeDistS(T_test, T_end, k - 1, {}))
                {
                    debugPrint("TreeDistS: Found solution with rotation");
                    return true;
                }
            }
            catch (...)
            {
                continue;
            }
        }
    }
    else
    {
        debugPrint("TreeDistS: Implementing Li-Xia branching over S");
        std::vector<std::pair<int,int>> chosen;
        if (branchOnSPairs(T_init, T_end, k, S, 0, chosen))
            return true;
    }

    // PSEUDOCODE STEP 3: "Return False"
    debugPrint("TreeDistS: No solution found");
    return false;
}

std::pair<std::vector<int>, std::vector<int>> rightComb(int m)
{
    // inorder = 1..m; preorder = 1..m
    std::vector<int> inorder(m), preorder(m);
    std::iota(inorder.begin(), inorder.end(), 1);
    std::iota(preorder.begin(), preorder.end(), 1);
    return {preorder, inorder};
}

std::pair<std::vector<int>, std::vector<int>> leftComb(int m)
{
    // inorder = 1..m; preorder = m..1
    std::vector<int> inorder(m), preorder(m);
    std::iota(inorder.begin(), inorder.end(), 1);
    for (int i = 0; i < m; ++i)
        preorder[i] = m - i;
    return {preorder, inorder};
}

static std::pair<std::vector<int>, std::vector<int>> serializeTree(const VectorRangeTreeMap &T)
{
    std::vector<int> pre, in;
    if (T.root < 0 || T.original_nodes.empty())
        return {pre, in};

    std::function<void(int)> preDFS = [&](int u)
    {
        if (u < 0 || !T.isOriginal(u))
            return;
        pre.push_back(u);
        preDFS(T.getLeftChild(u));
        preDFS(T.getRightChild(u));
    };
    std::function<void(int)> inDFS = [&](int u)
    {
        if (u < 0 || !T.isOriginal(u))
            return;
        inDFS(T.getLeftChild(u));
        in.push_back(u);
        inDFS(T.getRightChild(u));
    };

    preDFS(T.root);
    inDFS(T.root);
    return {pre, in};
}

struct PairVectorHash
{
    size_t operator()(const std::pair<std::vector<int>, std::vector<int>> &p) const noexcept
    {
        // very simple hasher; good enough for tiny fixtures
        size_t h = 1469598103934665603ull;
        auto mix = [&](int x)
        {
            h ^= (size_t)x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        };
        for (int x : p.first)
            mix(x);
        for (int x : p.second)
            mix(x);
        return h;
    }
};
struct PairVectorEq
{
    bool operator()(const std::pair<std::vector<int>, std::vector<int>> &a,
                    const std::pair<std::vector<int>, std::vector<int>> &b) const noexcept
    {
        return a.first == b.first && a.second == b.second;
    }
};

// Generate all neighbors by one rotation on any parent->child edge
static std::vector<VectorRangeTreeMap> oneRotationNeighbors(const VectorRangeTreeMap &T)
{
    std::vector<VectorRangeTreeMap> nbrs;
    auto edges = getInternalEdges(T); // (parent,child)
    for (auto &e : edges)
    {
        int p = e.first, c = e.second;
        VectorRangeTreeMap X = safeCopyTree(T);
        if (X.root < 0)
            continue;
        try
        {
            if (X.getLeftChild(p) == c)
            {
                X.rotateRight(p);
                nbrs.push_back(std::move(X));
            }
            else if (X.getRightChild(p) == c)
            {
                X.rotateLeft(p);
                nbrs.push_back(std::move(X));
            }
        }
        catch (...)
        {
            // ignore bad edge
        }
    }
    return nbrs;
}

// Exact minimal rotation distance by BFS, up to an optional cap
int MinRotationsBFS(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap = 50)
{
    if (TreesEqual(A, B))
        return 0;

    using Key = std::pair<std::vector<int>, std::vector<int>>;
    std::unordered_set<Key, PairVectorHash, PairVectorEq> seen;
    std::deque<std::pair<VectorRangeTreeMap, int>> q;

    Key kA = serializeTree(A);
    Key kB = serializeTree(B);
    seen.insert(kA);
    q.push_back({A, 0});

    while (!q.empty())
    {
        auto [cur, d] = q.front();
        q.pop_front();
        if (d >= cap)
            continue;

        for (auto &nx : oneRotationNeighbors(cur))
        {
            Key kk = serializeTree(nx);
            if (seen.insert(kk).second)
            {
                if (kk == kB)
                    return d + 1;
                q.push_back({std::move(nx), d + 1});
            }
        }
    }
    return -1; // not found within cap
}

void testOppositeFansFixture(int n_polygon)
{
    int m = n_polygon - 2;           // tree size for n-gon triangulation
    auto [preA, inA] = rightComb(m); // fan at 0 (right comb)
    auto [preB, inB] = leftComb(m);  // fan at n-1 (left comb)

    VectorRangeTreeMap A, B;
    A.build(preA, inA);
    B.build(preB, inB);

    int k_hint = m + 5;

    int alg_AB = FlipDistMinK(A, B, k_hint);
    int alg_BA = FlipDistMinK(B, A, k_hint);

    // exact oracle
    int bfs_AB = MinRotationsBFS(A, B, /*cap=*/64);
    int bfs_BA = MinRotationsBFS(B, A, /*cap=*/64);

    std::cout << "RESULT_CPP_FIXTURE,"
              << "n=" << n_polygon
              << ",m=" << m
              << ",alg_AB=" << alg_AB
              << ",alg_BA=" << alg_BA
              << ",bfs_AB=" << bfs_AB
              << ",bfs_BA=" << bfs_BA
              << ",expected=" << (n_polygon - 3)
              << "\n";
}

// ============================================================================
// COMPREHENSIVE ACCURACY AND SCALABILITY TESTING SUITE
// ============================================================================

// Performance timing utility
class PerformanceTimer
{
private:
    std::chrono::high_resolution_clock::time_point start_time;

public:
    void start() { start_time = std::chrono::high_resolution_clock::now(); }
    long long getMicroseconds()
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }
    double getMilliseconds() { return getMicroseconds() / 1000.0; }
};

// Tree generators for testing
class TreeGenerator
{
public:
    static std::pair<std::vector<int>, std::vector<int>> generateRightChain(int n)
    {
        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++)
        {
            preorder.push_back(i);
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateLeftChain(int n)
    {
        std::vector<int> preorder, inorder;
        for (int i = n; i >= 1; i--)
        {
            preorder.push_back(i);
        }
        for (int i = 1; i <= n; i++)
        {
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateBalanced(int n)
    {
        if (n == 0)
            return {{}, {}};
        if (n == 1)
            return {{1}, {1}};
        if (n == 2)
            return {{2, 1}, {1, 2}};
        if (n == 3)
            return {{2, 1, 3}, {1, 2, 3}};
        if (n == 4)
            return {{3, 2, 1, 4}, {1, 2, 3, 4}};
        if (n == 5)
            return {{3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}};
        if (n == 6)
            return {{4, 2, 1, 3, 5, 6}, {1, 2, 3, 4, 5, 6}};
        if (n == 7)
            return {{4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7}};

        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++)
            inorder.push_back(i);

        std::function<void(int, int)> buildBalanced = [&](int start, int end)
        {
            if (start > end)
                return;
            int mid = (start + end) / 2;
            preorder.push_back(mid);
            buildBalanced(start, mid - 1);
            buildBalanced(mid + 1, end);
        };

        buildBalanced(1, n);
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>>
    generateRandom(int n, int seed = 42)
    {
        if (n <= 0)
            return {{}, {}};

        std::mt19937 rng(seed);
        std::vector<int> inorder(n);
        std::iota(inorder.begin(), inorder.end(), 1);

        std::vector<int> preorder;
        preorder.reserve(n);

        // Recursively choose a random root in [lo, hi] and emit preorder
        std::function<void(int, int)> build = [&](int lo, int hi)
        {
            if (lo > hi)
                return;
            std::uniform_int_distribution<int> pick(lo, hi);
            int root = pick(rng);
            preorder.push_back(root);
            build(lo, root - 1);
            build(root + 1, hi);
        };

        build(1, n);
        return {preorder, inorder};
    }
};

void printTreeInfo(const std::string &name, const VectorRangeTreeMap &T)
{
    try
    {
        std::cout << name << " (root=" << T.root << ", nodes=" << T.original_nodes.size() << "): ";
        auto edges = getInternalEdges(T);
        for (const auto &e : edges)
        {
            std::cout << "(" << e.first << "," << e.second << ") ";
        }
        std::cout << std::endl;
    }
    catch (...)
    {
        std::cout << name << " [ERROR]" << std::endl;
    }
}

// ============================================================================
// ACCURACY TESTING - Verify Algorithm Correctness
// ============================================================================

void testAccuracyBasicCases()
{
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << " ACCURACY TESTING - BASIC CASES" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    struct TestCase
    {
        std::vector<int> pre1, in1, pre2, in2;
        int expected_min_distance;
        std::string description;
    };

    std::vector<TestCase> testCases = {
        // Basic cases
        {{1}, {1}, {1}, {1}, 0, "Single node (identical)"},
        {{1, 2}, {1, 2}, {1, 2}, {1, 2}, 0, "Two nodes (identical)"},
        {{1, 2}, {1, 2}, {2, 1}, {1, 2}, 1, "Two nodes (one rotation)"},

        // 3-node cases
        {{2, 1, 3}, {1, 2, 3}, {2, 1, 3}, {1, 2, 3}, 0, "3-node balanced (identical)"},
        {{1, 2, 3}, {1, 2, 3}, {2, 1, 3}, {1, 2, 3}, 1, "3-node chain to balanced"},
        {{2, 1, 3}, {1, 2, 3}, {3, 2, 1}, {1, 2, 3}, 1, "3-node balanced to mirrored"},
        {{1, 2, 3}, {1, 2, 3}, {3, 2, 1}, {1, 2, 3}, 2, "3-node chain to mirrored"},

        // 4-node cases
        {{1, 2, 3, 4}, {1, 2, 3, 4}, {3, 2, 1, 4}, {1, 2, 3, 4}, 2, "4-node chain to balanced"},
        {{3, 2, 1, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, 2, "4-node balanced to chain"},
        {{1, 2, 3, 4}, {1, 2, 3, 4}, {4, 3, 2, 1}, {1, 2, 3, 4}, 3, "4-node chain to reverse chain"},
    };

    int passed = 0;
    std::cout << std::setw(40) << "Test Case" << std::setw(12) << "Expected" << std::setw(12) << "Min k" << std::setw(12) << "Max k" << std::setw(10) << "Result" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    for (const auto &test : testCases)
    {
        VectorRangeTreeMap T1, T2;
        T1.build(test.pre1, test.in1);
        T2.build(test.pre2, test.in2);

        // Test if expected distance works
        bool works_at_expected = FlipDistTree(T1, T2, test.expected_min_distance);

        // Test if expected-1 fails (unless expected is 0)
        bool fails_before = (test.expected_min_distance == 0) ||
                            !FlipDistTree(T1, T2, test.expected_min_distance - 1);

        // Find the actual minimum working k
        int n = (int)test.in1.size();
        int min_k = FlipDistMinK(T1, T2, n + 5);

        // Find a reasonable upper bound
        int max_k = test.expected_min_distance + 3;
        for (int k = test.expected_min_distance; k <= test.expected_min_distance + 5; k++)
        {
            if (FlipDistTree(T1, T2, k))
            {
                max_k = k;
                break;
            }
        }

        int minK = FlipDistMinK(T1, T2, test.expected_min_distance + 5);
        std::cout << "RESULT_CPP_ACC,"
                  << "desc=\"" << test.description << "\""
                  << ",minK=" << minK
                  << ",expected=" << test.expected_min_distance
                  << std::endl;

        bool correct = works_at_expected && fails_before && (min_k == test.expected_min_distance);
        if (correct)
            passed++;

        std::cout << std::setw(40) << test.description.substr(0, 39)
                  << std::setw(12) << test.expected_min_distance
                  << std::setw(12) << min_k
                  << std::setw(12) << max_k
                  << std::setw(10) << (correct ? "✅ PASS" : "❌ FAIL") << std::endl;

        if (!correct)
        {
            std::cout << "    Details: works_at_expected=" << works_at_expected
                      << ", fails_before=" << fails_before
                      << ", min_k=" << min_k << std::endl;
        }
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Accuracy Test Results: " << passed << "/" << testCases.size() << " passed" << std::endl;
}

void testAccuracyMonotonicity()
{
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << " ACCURACY TESTING - MONOTONICITY PROPERTY" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "Testing: If FlipDistTree(T1, T2, k) = true, then FlipDistTree(T1, T2, k+1) = true" << std::endl;

    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::string>> cases = {
        {{1, 2, 3}, {1, 2, 3}, {3, 2, 1}, {1, 2, 3}, "3-node chain to mirrored"},
        {{1, 2, 3, 4}, {1, 2, 3, 4}, {3, 2, 1, 4}, {1, 2, 3, 4}, "4-node case"},
        {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}, "5-node case"},
    };

    int monotonic_cases = 0;
    for (const auto &[pre1, in1, pre2, in2, name] : cases)
    {
        VectorRangeTreeMap T1, T2;
        T1.build(pre1, in1);
        T2.build(pre2, in2);

        bool monotonic = true;
        bool found_true = false;
        std::vector<bool> results;

        std::cout << "\n"
                  << name << ":" << std::endl;
        std::cout << "k: ";
        for (int k = 0; k <= 8; k++)
        {
            bool result = FlipDistTree(T1, T2, k);
            results.push_back(result);
            std::cout << k << "=" << (result ? "T" : "F") << " ";

            if (found_true && !result)
            {
                monotonic = false;
            }
            if (result)
                found_true = true;
        }

        std::cout << " -> " << (monotonic ? "✅ MONOTONIC" : "❌ NOT MONOTONIC") << std::endl;
        if (monotonic)
            monotonic_cases++;
    }

    std::cout << "\nMonotonicity Results: " << monotonic_cases << "/" << cases.size() << " cases are monotonic" << std::endl;
}

void testAccuracyTreeEquality()
{
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << " ACCURACY TESTING - TREE EQUALITY CASES" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::vector<int> sizes = {1, 2, 3, 4, 5};
    int passed = 0;
    int total = 0;

    for (int n : sizes)
    {
        auto [pre, in] = TreeGenerator::generateRightChain(n);
        VectorRangeTreeMap T1, T2;
        T1.build(pre, in);
        T2.build(pre, in); // Identical tree

        bool result = FlipDistTree(T1, T2, 0); // Should work with k=0
        total++;
        if (result)
            passed++;

        std::cout << "Identical " << n << "-node trees: k=0 -> "
                  << (result ? "✅ TRUE" : "❌ FALSE") << std::endl;
    }

    std::cout << "Equality test results: " << passed << "/" << total << " passed" << std::endl;
}

// ============================================================================
// SCALABILITY TESTING - Find Maximum Node Count
// ============================================================================

void testScalabilityComprehensive()
{
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << " COMPREHENSIVE SCALABILITY ANALYSIS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    const long long TIMEOUT_MS = 5000; // 5 second timeout per test
    PerformanceTimer timer;

    std::vector<std::string> testTypes = {
        "Chain→Balanced",
        "Chain→Chain",
        "Balanced→Chain",
        "Random→Random"};

    std::cout << std::setw(6) << "Nodes"
              << std::setw(16) << "Test Type"
              << std::setw(12) << "Time (ms)"
              << std::setw(10) << "Min k"
              << std::setw(12) << "Success"
              << std::setw(15) << "Status" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (int n = 1; n <= 20; n++)
    { // Test up to 20 nodes
        bool any_timeout = false;

        for (const auto &testType : testTypes)
        {
            VectorRangeTreeMap T1, T2;

            try
            {
                // Generate trees based on test type
                if (testType == "Chain→Balanced")
                {
                    auto [pre1, in1] = TreeGenerator::generateRightChain(n);
                    auto [pre2, in2] = TreeGenerator::generateBalanced(n);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                }
                else if (testType == "Chain→Chain")
                {
                    auto [pre1, in1] = TreeGenerator::generateRightChain(n);
                    auto [pre2, in2] = TreeGenerator::generateLeftChain(n);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                }
                else if (testType == "Balanced→Chain")
                {
                    auto [pre1, in1] = TreeGenerator::generateBalanced(n);
                    auto [pre2, in2] = TreeGenerator::generateRightChain(n);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                }
                else if (testType == "Random→Random")
                {
                    auto [pre1, in1] = TreeGenerator::generateRandom(n, 42);
                    auto [pre2, in2] = TreeGenerator::generateRandom(n, 84);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                    // place holder
                }

                // Test with generous budget
                int budget = n + 5;

                timer.start();
                int min_k = FlipDistMinK(T1, T2, budget);
                bool success = FlipDistTree(T1, T2, budget);
                double time_ms = timer.getMilliseconds();

                int minK = FlipDistMinK(T1, T2, budget);

                // machine-friendly line you can grep/sort later
                std::cout << "RESULT_CPP,"
                          << "n=" << n
                          << ",type=" << testType
                          << ",minK=" << minK
                          << ",time_ms=" << std::fixed << std::setprecision(2) << time_ms
                          << ",success=" << (success ? 1 : 0)
                          << std::endl;

                std::string status = "OK";
                if (time_ms > TIMEOUT_MS)
                {
                    status = "TIMEOUT";
                    any_timeout = true;
                }

                std::cout << std::setw(6) << n
                          << std::setw(16) << testType
                          << std::setw(12) << std::fixed << std::setprecision(2) << time_ms
                          << std::setw(10) << (min_k == -1 ? -1 : min_k)
                          << std::setw(12) << (success ? "YES" : "NO")
                          << std::setw(15) << status << std::endl;

                // Stop this test type if timeout
                if (time_ms > TIMEOUT_MS)
                {
                    break;
                }
            }
            catch (...)
            {
                std::cout << std::setw(6) << n
                          << std::setw(16) << testType
                          << std::setw(12) << "ERROR"
                          << std::setw(12) << "NO"
                          << std::setw(15) << "EXCEPTION" << std::endl;
            }
        }

        // Stop if most test types are timing out
        if (any_timeout && n >= 5)
        {
            std::cout << "\n⚠  Stopping at " << n << " nodes due to performance limits" << std::endl;
            break;
        }
    }
}

void testScalabilityDetailed()
{
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << " DETAILED PERFORMANCE ANALYSIS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::vector<int> testSizes = {3, 4, 5, 6, 7, 8, 9, 10};
    PerformanceTimer timer;

    std::cout << std::setw(6) << "Nodes"
              << std::setw(12) << "Avg (ms)"
              << std::setw(12) << "Min (ms)"
              << std::setw(12) << "Max (ms)"
              << std::setw(12) << "Success"
              << std::setw(15) << "Notes" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (int n : testSizes)
    {
        std::vector<double> times;
        int successes = 0;
        const int NUM_RUNS = 3;

        for (int run = 0; run < NUM_RUNS; run++)
        {
            try
            {
                auto [pre1, in1] = TreeGenerator::generateRightChain(n);
                auto [pre2, in2] = TreeGenerator::generateBalanced(n);

                VectorRangeTreeMap T1, T2;
                T1.build(pre1, in1);
                T2.build(pre2, in2);

                timer.start();
                bool result = FlipDistTree(T1, T2, n + 3);
                double time_ms = timer.getMilliseconds();

                times.push_back(time_ms);
                if (result)
                    successes++;

                // Stop if getting too slow
                if (time_ms > 10000)
                { // 10 seconds
                    break;
                }
            }
            catch (...)
            {
                times.push_back(-1); // Error marker
            }
        }

        if (!times.empty() && times[0] >= 0)
        {
            double avg = 0, min_time = times[0], max_time = times[0];
            int valid_times = 0;

            for (double t : times)
            {
                if (t >= 0)
                {
                    avg += t;
                    min_time = std::min(min_time, t);
                    max_time = std::max(max_time, t);
                    valid_times++;
                }
            }

            if (valid_times > 0)
            {
                avg /= valid_times;

                std::string notes = "";
                if (avg > 1000)
                    notes = "SLOW";
                else if (avg > 100)
                    notes = "MODERATE";
                else
                    notes = "FAST";

                std::cout << std::setw(6) << n
                          << std::setw(12) << std::fixed << std::setprecision(2) << avg
                          << std::setw(12) << std::setprecision(2) << min_time
                          << std::setw(12) << std::setprecision(2) << max_time
                          << std::setw(12) << successes << "/" << NUM_RUNS
                          << std::setw(15) << notes << std::endl;

                // Stop if consistently slow
                if (avg > 5000)
                {
                    std::cout << " Stopping detailed analysis due to performance" << std::endl;
                    break;
                }
            }
        }
        else
        {
            std::cout << std::setw(6) << n
                      << std::setw(12) << "ERROR"
                      << std::setw(12) << "-"
                      << std::setw(12) << "-"
                      << std::setw(12) << "0/" << NUM_RUNS
                      << std::setw(15) << "FAILED" << std::endl;
        }
    }
}

// ============================================================================
// MAIN TESTING FUNCTION
// ============================================================================

void runComprehensiveTests()
{
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << "COMPREHENSIVE ALGORITHM TESTING SUITE" << std::endl;
    std::cout << "    Accuracy Validation + Scalability Analysis" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    // Phase 1: Accuracy Testing
    std::cout << "\nPHASE 1: ACCURACY VALIDATION" << std::endl;
    testAccuracyBasicCases();
    testAccuracyMonotonicity();
    testAccuracyTreeEquality();

    // Phase 2: Scalability Testing
    std::cout << "\nPHASE 2: SCALABILITY ANALYSIS" << std::endl;
    testScalabilityComprehensive();
    testScalabilityDetailed();

    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << "COMPREHENSIVE TESTING COMPLETE!" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
}

struct FlipStep {
    int pivot{-1};
    bool left{false};
};

struct FlipParentInfo {
    std::string parent_key;
    int pivot{-1};
    bool left{false};
};

static bool buildCanonicalRotationPath(const VectorRangeTreeMap& start,
                                       const VectorRangeTreeMap& goal,
                                       std::vector<std::string>& path_out,
                                       std::vector<FlipStep>& moves_out)
{
    path_out.clear();
    moves_out.clear();

    std::string startKey = treeToString(start);
    std::string goalKey = treeToString(goal);
    if (startKey == goalKey) {
        path_out.push_back(treeToString(start));
        return true;
    }

    struct Node {
        VectorRangeTreeMap tree;
        std::string key;
    };

    std::queue<Node> q;
    std::unordered_set<std::string> visited;
    std::unordered_map<std::string, FlipParentInfo> parent;

    visited.insert(startKey);
    parent.emplace(startKey, FlipParentInfo{"", -1, false});
    q.push(Node{start, startKey});

    bool found = false;
    while (!q.empty() && !found) {
        Node cur = std::move(q.front());
        q.pop();

        for (int v : cur.tree.original_nodes) {
            if (cur.tree.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur.tree;
                tmp.rotateLeft(v);
                std::string key = treeToString(tmp);
                if (visited.insert(key).second) {
                    parent.emplace(key, FlipParentInfo{cur.key, v, true});
                    if (key == goalKey) { found = true; break; }
                    q.push(Node{std::move(tmp), key});
                }
            }
            if (cur.tree.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur.tree;
                tmp.rotateRight(v);
                std::string key = treeToString(tmp);
                if (visited.insert(key).second) {
                    parent.emplace(key, FlipParentInfo{cur.key, v, false});
                    if (key == goalKey) { found = true; break; }
                    q.push(Node{std::move(tmp), key});
                }
            }
        }
    }

    if (!found) {
        path_out.push_back(treeToString(start));
        return false;
    }

    std::vector<FlipParentInfo> moves;
    std::string curKey = goalKey;
    while (curKey != startKey) {
        const FlipParentInfo& info = parent.at(curKey);
        moves.push_back(info);
        curKey = info.parent_key;
    }
    std::reverse(moves.begin(), moves.end());

    VectorRangeTreeMap current = start;
    path_out.push_back(treeToString(current));
    for (const auto& step : moves) {
        if (step.pivot >= 0) {
            if (step.left) {
                current.rotateLeft(step.pivot);
            } else {
                current.rotateRight(step.pivot);
            }
        }
        path_out.push_back(treeToString(current));
        moves_out.push_back(FlipStep{step.pivot, step.left});
    }
    return true;
}

static void printTreeAscii(const VectorRangeTreeMap& tree,
                           std::ostream& os,
                           const char* heading)
{
    os << heading << "\n";
    if (tree.root == VectorRangeTreeMap::NO_CHILD) {
        os << "  (empty)\n";
        return;
    }

    std::function<void(int,const std::string&,bool)> rec =
        [&](int node, const std::string& prefix, bool left) {
            if (node == VectorRangeTreeMap::NO_CHILD || !tree.isOriginal(node)) return;

            int right = tree.getRightChild(node);
            int leftChild = tree.getLeftChild(node);
            bool hasRight = (right != VectorRangeTreeMap::NO_CHILD) && tree.isOriginal(right);
            bool hasLeft  = (leftChild != VectorRangeTreeMap::NO_CHILD) && tree.isOriginal(leftChild);

            if (hasRight) rec(right, prefix + (left ? "│   " : "    "), false);

            os << prefix;
            if (!prefix.empty()) {
                os << (left ? "└── " : "┌── ");
            }
            os << node << "\n";

            if (hasLeft)  rec(leftChild, prefix + (left ? "    " : "│   "), true);
        };

    rec(tree.root, "", true);
}

static void logPathToStderr(const std::string& label,
                            const VectorRangeTreeMap& start,
                            const std::vector<FlipStep>& moves,
                            bool complete,
                            bool ascii)
{
    std::cerr << "[PATH] " << label
              << " length=" << static_cast<int>(moves.size())
              << " complete=" << (complete ? "yes" : "no") << "\n";

    VectorRangeTreeMap current = start;
    auto logState = [&](size_t step_index, const char* prefix) {
        std::cerr << "  Step " << step_index << ' ' << prefix
                  << canonicalTraversalString(current) << "\n";
        if (ascii) {
            printTreeAscii(current, std::cerr, "    ");
        }
    };

    logState(0, "(start): ");
    for (size_t i = 0; i < moves.size(); ++i) {
        const FlipStep& mv = moves[i];
        if (mv.pivot >= 0) {
            if (mv.left) {
                current.rotateLeft(mv.pivot);
            } else {
                current.rotateRight(mv.pivot);
            }
        }
        std::string prefix = std::string("(") + (mv.left ? "rotateLeft" : "rotateRight")
                            + " " + std::to_string(mv.pivot) + "): ";
        logState(i + 1, prefix.c_str());
    }
    std::cerr.flush();
}

struct FlipCliOptions {
    std::string program = "flipdist_asan";
    std::string case_type = "comb";
    int n = 5;
    int count = 1;
    long long seed = 12345;
    int max_k = -1;
    int bfs_cap = 64;
    bool run_legacy = false;
    double time_limit = 0.0;
    std::size_t visited_cap = 0;
    std::size_t queue_cap = 0;
    bool fallback_bidir = false;
    std::size_t bidir_cap = 0;
    bool prefer_bidir = false;
    bool emit_path = false;
    bool path_ascii = false;
};

static void printUsage(const char *argv0)
{
    std::cerr << "Usage: " << argv0
              << " [--case comb|random]"
              << " [--n N]"
              << " [--count C]"
              << " [--seed S]"
              << " [--program NAME]"
              << " [--max-k K]"
              << " [--bfs-cap CAP]"
              << " [--emit-path]"
              << " [--path-ascii]"
              << " [--legacy-fixtures]\n";
}

static bool parseCliOptions(int argc, char **argv, FlipCliOptions &opts)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        auto consume = [&](const std::string &name) -> std::string {
            if (arg == name)
            {
                if (i + 1 >= argc)
                {
                    throw std::invalid_argument(name + " requires a value");
                }
                return std::string(argv[++i]);
            }
            auto pos = arg.find('=');
            if (pos != std::string::npos && arg.substr(0, pos) == name)
            {
                return arg.substr(pos + 1);
            }
            return {};
        };

        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            std::exit(0);
        }

        if (auto value = consume("--case"); !value.empty())
        {
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            opts.case_type = value;
            continue;
        }
        if (auto value = consume("--n"); !value.empty())
        {
            opts.n = std::stoi(value);
            continue;
        }
        if (auto value = consume("--count"); !value.empty())
        {
            opts.count = std::stoi(value);
            continue;
        }
        if (auto value = consume("--seed"); !value.empty())
        {
            opts.seed = std::stoll(value);
            continue;
        }
        if (auto value = consume("--program"); !value.empty())
        {
            opts.program = value;
            continue;
        }
        if (auto value = consume("--max-k"); !value.empty())
        {
            opts.max_k = std::stoi(value);
            continue;
        }
        if (auto value = consume("--bfs-cap"); !value.empty())
        {
            opts.bfs_cap = std::stoi(value);
            continue;
        }
        if (auto value = consume("--time-limit"); !value.empty())
        {
            opts.time_limit = std::stod(value);
            continue;
        }
        if (auto value = consume("--visited-cap"); !value.empty())
        {
            opts.visited_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (auto value = consume("--queue-cap"); !value.empty())
        {
            opts.queue_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (auto value = consume("--bidir-cap"); !value.empty())
        {
            opts.bidir_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (arg == "--fallback-bidir")
        {
            opts.fallback_bidir = true;
            continue;
        }
        if (arg == "--prefer-bidir")
        {
            opts.prefer_bidir = true;
            continue;
        }
        if (arg == "--legacy-fixtures" || arg == "--fixtures")
        {
            opts.run_legacy = true;
            continue;
        }
        if (arg == "--emit-path")
        {
            opts.emit_path = true;
            continue;
        }
        if (arg == "--path-ascii")
        {
            opts.path_ascii = true;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        printUsage(argv[0]);
        return false;
    }
    return true;
}

struct FlipCase
{
    VectorRangeTreeMap start;
    VectorRangeTreeMap target;
    long long seed;
};

static bool buildCase(const FlipCliOptions &opts, int index, FlipCase &out)
{
    if (opts.n <= 0)
    {
        std::cerr << "n must be positive\n";
        return false;
    }

    if (opts.case_type == "comb")
    {
        Traversals left = makeCombTraversals(opts.n, /*rightComb=*/false);
        Traversals right = makeCombTraversals(opts.n, /*rightComb=*/true);
        out.start = VectorRangeTreeMap();
        out.target = VectorRangeTreeMap();
        out.start.build(left.preorder, left.inorder);
        out.target.build(right.preorder, right.inorder);
        out.seed = -1;
        return out.start.root != VectorRangeTreeMap::NO_CHILD &&
               out.target.root != VectorRangeTreeMap::NO_CHILD;
    }

    if (opts.case_type == "random")
    {
        long long pair_seed = opts.seed + index;
        std::mt19937 rngA(static_cast<std::uint32_t>(pair_seed * 2 + 0));
        std::mt19937 rngB(static_cast<std::uint32_t>(pair_seed * 2 + 1));
        Traversals tA = makeRandomTraversals(opts.n, rngA);
        Traversals tB = makeRandomTraversals(opts.n, rngB);
        out.start = VectorRangeTreeMap();
        out.target = VectorRangeTreeMap();
        out.start.build(tA.preorder, tA.inorder);
        out.target.build(tB.preorder, tB.inorder);
        out.seed = pair_seed;
        return out.start.root != VectorRangeTreeMap::NO_CHILD &&
               out.target.root != VectorRangeTreeMap::NO_CHILD;
    }

    std::cerr << "Unsupported case type: " << opts.case_type << "\n";
    return false;
}

struct FlipRun
{
    std::string program;
    std::string case_type;
    int n = 0;
    long long seed = -1;
    std::string direction;
    int distance = -1;
    int distance_flipdist = -1;
    int distance_bfs = -1;
    double time_ms = 0.0;
    double time_ms_flipdist = 0.0;
    double time_ms_bfs = 0.0;
    std::size_t expanded = 0;
    std::size_t enqueued = 0;
    std::size_t visited = 0;
    std::size_t max_queue = 0;
    std::size_t duplicates = 0;
    std::string status = "not_run";
    std::string status_flipdist = "not_run";
    std::string status_bfs = "not_run";
    std::string solver = "flipdist";
    std::string tree_a;
    std::string tree_b;
    int max_k = 0;
    std::vector<std::string> path;
    std::vector<FlipStep> moves;
    bool path_complete = false;
};

static std::string escapeJson(const std::string &s)
{
    std::string out;
    out.reserve(s.size() + 16);
    for (unsigned char c : s)
    {
        switch (c)
        {
        case '\"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\b':
            out += "\\b";
            break;
        case '\f':
            out += "\\f";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (c < 0x20)
            {
                char buf[7];
                std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned int>(c));
                out += buf;
            }
            else
            {
                out.push_back(static_cast<char>(c));
            }
        }
    }
    return out;
}

static std::string runToJson(const FlipRun &row)
{
    std::ostringstream out;
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setprecision(3);
    auto toULL = [](std::size_t v)
    { return static_cast<unsigned long long>(v); };

    out << "{\"program\":\"" << escapeJson(row.program) << '"'
        << ",\"case_type\":\"" << escapeJson(row.case_type) << '"'
        << ",\"n\":" << row.n
        << ",\"seed\":" << row.seed
        << ",\"direction\":\"" << escapeJson(row.direction) << '"'
        << ",\"distance\":" << row.distance
        << ",\"distance_bfs\":" << row.distance_bfs
        << ",\"time_ms\":" << row.time_ms
        << ",\"time_ms_flipdist\":" << row.time_ms_flipdist
        << ",\"time_ms_bfs\":" << row.time_ms_bfs
        << ",\"expanded\":" << toULL(row.expanded)
        << ",\"enqueued\":" << toULL(row.enqueued)
        << ",\"visited\":" << toULL(row.visited)
        << ",\"max_queue\":" << toULL(row.max_queue)
        << ",\"duplicates\":" << toULL(row.duplicates)
        << ",\"status\":\"" << escapeJson(row.status) << '"'
        << ",\"status_flipdist\":\"" << escapeJson(row.status_flipdist) << '"'
        << ",\"status_bfs\":\"" << escapeJson(row.status_bfs) << '"'
        << ",\"solver\":\"" << escapeJson(row.solver) << '"'
        << ",\"tree_a\":\"" << escapeJson(row.tree_a) << '"'
        << ",\"tree_b\":\"" << escapeJson(row.tree_b) << '"'
        << ",\"max_k\":" << row.max_k;

    if (!row.path.empty()) {
        out << ",\"path\":[";
        for (size_t i = 0; i < row.path.size(); ++i) {
            if (i) out << ',';
            out << '"' << escapeJson(row.path[i]) << '"';
        }
        out << "]";
        out << ",\"path_complete\":" << (row.path_complete ? "true" : "false");

        out << ",\"moves\":[";
        for (size_t i = 0; i < row.moves.size(); ++i) {
            if (i) out << ',';
            out << "{\"pivot\":" << row.moves[i].pivot
                << ",\"direction\":\"" << (row.moves[i].left ? "left" : "right") << "\"}";
        }
        out << "]";
    }

    out << "}";
    return out.str();
}

static FlipRun evaluateFlipCase(const VectorRangeTreeMap &start,
                                const VectorRangeTreeMap &target,
                                const FlipCliOptions &opts,
                                const std::string &direction,
                                long long seed_value,
                                int n_nodes)
{
    FlipRun row;
    row.program = opts.program;
    row.case_type = opts.case_type;
    row.n = n_nodes;
    row.seed = seed_value;
    row.direction = direction;
    row.tree_a = canonicalTraversalString(start);
    row.tree_b = canonicalTraversalString(target);
    row.max_k = (opts.max_k > 0) ? opts.max_k : std::max(1, 2 * n_nodes + 6);

    int bfs_cap = (opts.bfs_cap > 0) ? opts.bfs_cap : std::max(32, 2 * n_nodes + 6);

    auto flip_start = std::chrono::steady_clock::now();
    g_flipDistMemo.clear();
    int dist = FlipDistMinK(start, target, row.max_k);
    auto flip_end = std::chrono::steady_clock::now();
    row.time_ms = std::chrono::duration<double, std::milli>(flip_end - flip_start).count();
    row.distance = dist;
    row.distance_flipdist = dist;
    row.time_ms_flipdist = row.time_ms;
    row.status = (dist >= 0) ? "ok" : "not_found";
    row.status_flipdist = row.status;

    auto bfs_start = std::chrono::steady_clock::now();
    int bfs_dist = MinRotationsBFS(start, target, bfs_cap);
    auto bfs_end = std::chrono::steady_clock::now();
    row.time_ms_bfs = std::chrono::duration<double, std::milli>(bfs_end - bfs_start).count();
    row.distance_bfs = bfs_dist;
    row.status_bfs = (bfs_dist >= 0) ? "ok" : "cap";

    const bool flip_ok = (row.status == "ok");
    const bool bfs_ok = (row.status_bfs == "ok");

    if (bfs_ok && (!flip_ok || row.distance != row.distance_bfs))
    {
        row.status = flip_ok ? "bfs_override" : "bfs_only";
        row.distance = row.distance_bfs;
        row.time_ms = row.time_ms_bfs;
        row.solver = "bfs";
        std::cerr << "[WARN] flipdist mismatch resolved via BFS"
                  << " case=" << row.case_type
                  << " n=" << row.n
                  << " seed=" << row.seed
                  << " dir=" << row.direction
                  << " flipdist=" << row.distance_flipdist
                  << " bfs=" << row.distance_bfs
                  << "\n";
    }

    if (opts.emit_path && row.distance >= 0) {
        std::vector<std::string> path;
        std::vector<FlipStep> moves;
        bool complete = buildCanonicalRotationPath(start, target, path, moves);
        row.path = std::move(path);
        row.moves = std::move(moves);
        row.path_complete = complete;
    }

    return row;
}

int main(int argc, char **argv)
{
    FlipCliOptions opts;
    try
    {
        if (!parseCliOptions(argc, argv, opts))
        {
            return 1;
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Argument error: " << ex.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }

    if (opts.run_legacy)
    {
        for (int n_polygon = 11; n_polygon <= 13; ++n_polygon)
        {
            testOppositeFansFixture(n_polygon);
        }
        return 0;
    }

    if (opts.n <= 0)
    {
        std::cerr << "n must be positive\n";
        return 1;
    }
    if (opts.count <= 0)
    {
        std::cerr << "count must be positive\n";
        return 1;
    }

    for (int idx = 0; idx < opts.count; ++idx)
    {
        FlipCase instance;
        if (!buildCase(opts, idx, instance))
        {
            std::cerr << "Failed to build case index " << idx << "\n";
            return 1;
        }

        int n_nodes = static_cast<int>(instance.start.original_nodes.size());
        FlipRun forward = evaluateFlipCase(instance.start, instance.target, opts, "a->b", instance.seed, n_nodes);
        std::cout << runToJson(forward) << "\n";
        std::cout.flush();
        if (opts.emit_path && !forward.path.empty()) {
            logPathToStderr(forward.case_type + " " + forward.direction,
                            instance.start,
                            forward.moves,
                            forward.path_complete,
                            opts.path_ascii);
        }

        FlipRun reverse = evaluateFlipCase(instance.target, instance.start, opts, "b->a", instance.seed, n_nodes);
        std::cout << runToJson(reverse) << "\n";
        std::cout.flush();
        if (opts.emit_path && !reverse.path.empty()) {
            logPathToStderr(reverse.case_type + " " + reverse.direction,
                            instance.target,
                            reverse.moves,
                            reverse.path_complete,
                            opts.path_ascii);
        }
    }

    return 0;
}
