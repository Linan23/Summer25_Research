#include "helpers.h"
#include "algorithm.h"
#include "memoization.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <functional>

// Helper Functions

// Definition: Collect internal parent->child edges over original nodes
// Parameters: T: tree to scan
// Returns: vector of (parent, child) edges
// Errors: returns empty on invalid trees or if traversal fails
std::vector<std::pair<int,int>> getInternalEdges(const VectorRangeTreeMap& T) {
    std::vector<std::pair<int,int>> edges;
    try {
        if (T.original_nodes.empty() || T.root < 0) {
            return edges;
        }

        std::function<void(int)> dfs = [&](int node) {
            if (node < 0 || !T.isOriginal(node)) return;
            try {
                int left = T.getLeftChild(node);
                int right = T.getRightChild(node);
                if (left >= 0 && T.isOriginal(left)) {
                    edges.emplace_back(node, left);
                    dfs(left);
                }
                if (right >= 0 && T.isOriginal(right)) {
                    edges.emplace_back(node, right);
                    dfs(right);
                }
            } catch (...) {
                return;
            }
        };

        if (T.isOriginal(T.root)) {
            dfs(T.root);
        }
    } catch (...) {
        edges.clear();
    }
    return edges;
}

// Definition: Count internal parent->child edges in a tree
// Parameters: T: tree to scan
// Returns: number of edges
// Errors: returns 0 on invalid trees or traversal failure
int countInternalEdges(const VectorRangeTreeMap& T) {
    if (T.original_nodes.empty()) {
        return 0;
    }
    return static_cast<int>(T.original_nodes.size()) - 1;
}

// Definition: Count edges in T_init that are not present in T_final (by range)
// Parameters: T_init: source tree; T_final: target tree
// Returns: conflict count
// Errors: returns 0 on invalid trees
int countConflictEdges(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final) {
    if (T_init.original_nodes.empty() || T_final.original_nodes.empty()) return 0;
    auto targetEdges = buildTargetSet(T_final);
    int conflicts = 0;
    for (int v : T_init.original_nodes) {
        if (!T_init.isOriginal(v)) continue;
        auto pr = T_init.getRange(v);
        for (int c : {T_init.getLeftChild(v), T_init.getRightChild(v)}) {
            if (c >= 0 && T_init.isOriginal(c)) {
                auto cr = T_init.getRange(c);
                RP e{pr.first, pr.second, cr.first, cr.second};
                if (!targetEdges.count(e)) {
                    conflicts++;
                }
            }
        }
    }
    return conflicts;
}

// Definition: Check if two edges share a node
// Parameters: e1: first edge; e2: second edge
// Returns: true if they share an endpoint, else false
// Errors: none
bool areAdjacent(const std::pair<int,int>& e1, const std::pair<int,int>& e2) {
    return e1.first == e2.first || e1.first == e2.second ||
           e1.second == e2.first || e1.second == e2.second;
}

struct EndpointEntry {
    int other = -1;
    std::pair<int,int> edge{-1,-1};
    bool boundary = false;
};

// Definition: Compute vertex count (leaf count) for endpoint indexing
// Parameters: T: tree to inspect
// Returns: number of polygon vertices (leaf count)
// Errors: returns 0 on invalid tree
static int vertexCount(const VectorRangeTreeMap& T) {
    if (T.root < 0 || !T.isOriginal(T.root)) return 0;
    auto r = T.getRange(T.root);
    return r.second + 1;
}

// Definition: Build endpoint index for diagonals plus boundary edges
// Parameters: T: tree; endpointMap/diagMap: output maps
// Returns: nothing
// Errors: outputs may be empty on invalid trees
static void buildEndpointIndex(const VectorRangeTreeMap& T,
                               std::unordered_map<int, std::vector<EndpointEntry>>& endpointMap) {
    endpointMap.clear();

    int vcount = vertexCount(T);
    if (vcount <= 0) return;

    auto edges = getInternalEdges(T);
    for (const auto& edge : edges) {
        int parent = edge.first;
        int child = edge.second;
        if (!T.isOriginal(child) || !T.isOriginal(parent)) continue;
        auto cr = T.getRange(child);
        int L = cr.first;
        int R = cr.second;
        endpointMap[L].push_back({R, edge, false});
        endpointMap[R].push_back({L, edge, false});
    }

    // Add boundary edges (cycle)
    for (int v = 0; v < vcount; v++) {
        int prev = (v - 1 + vcount) % vcount;
        int next = (v + 1) % vcount;
        endpointMap[v].push_back({prev, {-1,-1}, true});
        endpointMap[v].push_back({next, {-1,-1}, true});
    }
}

// Definition: Pick the two neighboring edges around an endpoint for a diagonal
// Parameters: endpoint: fixed endpoint; other: other endpoint of diagonal; vcount: vertex count; endpointMap: incident edges
// Returns: pair of edges (boundary edges use {-1,-1})
// Errors: returns boundary pair on missing data
static std::pair<std::pair<int,int>, std::pair<int,int>> pickNeighborPair(
        int endpoint,
        int other,
        int vcount,
        const std::unordered_map<int, std::vector<EndpointEntry>>& endpointMap) {

    auto it = endpointMap.find(endpoint);
    if (it == endpointMap.end() || it->second.empty() || vcount <= 0) {
        return {{-1,-1},{-1,-1}};
    }

    struct Item {
        int other;
        std::pair<int,int> edge;
        bool boundary;
        int dist;
    };

    std::vector<Item> items;
    items.reserve(it->second.size());
    for (const auto& entry : it->second) {
        int dist = (entry.other - endpoint + vcount) % vcount;
        items.push_back({entry.other, entry.edge, entry.boundary, dist});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.dist < b.dist;
    });

    int idx = -1;
    for (size_t i = 0; i < items.size(); i++) {
        if (!items[i].boundary && items[i].other == other) {
            idx = static_cast<int>(i);
            break;
        }
    }

    if (idx < 0) {
        return {{-1,-1},{-1,-1}};
    }

    int prev = (idx - 1 + (int)items.size()) % (int)items.size();
    int next = (idx + 1) % (int)items.size();

    auto e1 = items[prev].boundary ? std::make_pair(-1,-1) : items[prev].edge;
    auto e2 = items[next].boundary ? std::make_pair(-1,-1) : items[next].edge;
    return {e1, e2};
}

// Definition: Build a deep copy while refreshing preorder metadata
// Parameters: T: tree to copy
// Returns: copied tree, or an empty tree on failure
// Errors: returns empty on invalid trees or reconstruction failure
VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap& T) {
    try {
        if (T.root < 0 || !T.isOriginal(T.root) || T.original_nodes.empty()) {
            return VectorRangeTreeMap{};
        }

        VectorRangeTreeMap copy;
        copy.ranges = T.ranges;
        copy.edges = T.edges;
        copy.parents = T.parents;
        copy.root = T.root;
        copy.max_node_value = T.max_node_value;
        copy.original_inorder = T.original_inorder;
        copy.position_in_inorder_fast = T.position_in_inorder_fast;
        copy.original_nodes = T.original_nodes;
        copy.original_mask = T.original_mask;
        copy.fingerprint_valid = T.fingerprint_valid;
        copy.fingerprint_h1 = T.fingerprint_h1;
        copy.fingerprint_h2 = T.fingerprint_h2;
        copy.original_preorder.clear();
        copy.original_preorder.reserve(T.original_nodes.size());

        std::function<void(int)> buildPreorder = [&](int node) {
            if (node < 0 || !T.isOriginal(node)) return;
            copy.original_preorder.push_back(node);
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left)) buildPreorder(left);
            if (right >= 0 && T.isOriginal(right)) buildPreorder(right);
        };

        buildPreorder(T.root);
        if (copy.original_preorder.size() != T.original_nodes.size()) {
            return VectorRangeTreeMap{};
        }
        return copy;
    } catch (...) {
        return VectorRangeTreeMap{};
    }
}

// Definition: Check whether parent->child is a direct edge
// Parameters: T: tree to check; parent: parent node id; child: child node id
// Returns: true if parent has child as a direct left/right child
// Errors: returns false on invalid nodes or tree access failures
bool hasParentChildEdge(const VectorRangeTreeMap& T, int parent, int child) {
    try {
        if (!T.isOriginal(parent) || !T.isOriginal(child)) return false;
        int left = T.getLeftChild(parent);
        int right = T.getRightChild(parent);
        return (left == child) || (right == child);
    } catch (...) {
        return false;
    }
}

// Definition: Resolve an undirected edge to a valid parent->child orientation
// Parameters: T: tree to check; edge: undirected node pair
// Returns: oriented {parent, child} if edge exists, else {-1, -1}
// Errors: returns {-1, -1} on invalid nodes or if edge not present
std::pair<int,int> orientEdge(const VectorRangeTreeMap& T, const std::pair<int,int>& edge) {
    if (hasParentChildEdge(T, edge.first, edge.second)) {
        return edge;
    }
    if (hasParentChildEdge(T, edge.second, edge.first)) {
        return {edge.second, edge.first};
    }
    return {-1, -1};
}

// Definition: Find a rotatable edge that would insert a target edge
// Parameters: T_init: current tree; T_final: target tree
// Returns: {true, edge} if a free edge exists, otherwise {false, {-1,-1}}
// Errors: returns false on invalid trees or if rotations fail
std::pair<bool, std::pair<int,int>> findFreeEdge(const VectorRangeTreeMap& T_init,
                                                 const VectorRangeTreeMap& T_final) {
    int pivot = -1;
    bool leftRotation = false;
    if (!hasFreeEdge(T_init, T_final, pivot, leftRotation)) {
        return {false, {-1, -1}};
    }
    int child = leftRotation ? T_init.getRightChild(pivot) : T_init.getLeftChild(pivot);
    if (child < 0 || !T_init.isOriginal(child)) {
        return {false, {-1, -1}};
    }
    return {true, {pivot, child}};
}

// Definition: Try to decompose along a common edge shared by both trees
// Parameters: T_init: source tree; T_final: target tree; k: budget; handled: set true if a valid common edge is processed
// Returns: true if any decomposition succeeds within k
// Errors: returns false if no valid common edge is found or no split succeeds
bool tryCommonEdgeDecomposition(const VectorRangeTreeMap& T_init,
                                const VectorRangeTreeMap& T_final,
                                int k,
                                bool &handled) {
    handled = false;
    if (T_init.original_nodes.empty() || T_final.original_nodes.empty()) {
        return false;
    }

    auto targetEdges = buildTargetSet(T_final);

    const auto boundsMinSuccess = [](const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                    const Key128& key) {
        auto it = bounds_map.find(key);
        if (it == bounds_map.end()) return -1;
        return it->second.min_success;
    };
    const auto boundsMaxFail = [](const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                 const Key128& key) {
        auto it = bounds_map.find(key);
        if (it == bounds_map.end()) return -1;
        return it->second.max_fail;
    };

    std::vector<RP> commonEdges;
    for (int v : T_init.original_nodes) {
        if (!T_init.isOriginal(v)) continue;
        auto pr = T_init.getRange(v);
        for (int c : {T_init.getLeftChild(v), T_init.getRightChild(v)}) {
            if (c >= 0 && T_init.isOriginal(c)) {
                auto cr = T_init.getRange(c);
                RP e{pr.first, pr.second, cr.first, cr.second};
                if (targetEdges.count(e)) {
                    commonEdges.push_back(e);
                }
            }
        }
    }

    for (const auto& edge : commonEdges) {
        std::pair<int,int> parent_range{edge.ps, edge.pe};
        std::pair<int,int> child_range{edge.cs, edge.ce};

        auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(T_init, parent_range, child_range);
        auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(T_final, parent_range, child_range);

        if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes) {
            continue;
        }

        handled = true;
        int conf1 = countConflictEdges(A1, B1);
        int conf2 = countConflictEdges(A2, B2);
        if (conf1 + conf2 > k) {
            continue;
        }

        const Key128 p1_key = makeKeyPair(A1, B1);
        const Key128 p2_key = makeKeyPair(A2, B2);
        // Note: outer boundsMinSuccess is intentionally reused here.
        int lb1 = std::max(conf1, requiredBudgetFromBounds(g_kbounds, p1_key));
        int lb2 = std::max(conf2, requiredBudgetFromBounds(g_kbounds, p2_key));

        int min_k1 = std::max(0, lb1);
        int max_k1 = std::min(k, k - lb2);
        if (min_k1 > max_k1) {
            continue;
        }

        const int side1_min_success = boundsMinSuccess(g_kbounds, p1_key);
        const int side2_min_success = boundsMinSuccess(g_kbounds, p2_key);
        if (side1_min_success >= 0 && side2_min_success >= 0) {
            const int guaranteed_k1_min = std::max(min_k1, side1_min_success);
            const int guaranteed_k1_max = std::min(max_k1, k - side2_min_success);
            if (guaranteed_k1_min <= guaranteed_k1_max) {
                return true;
            }
        }

        int side1_fail_until = std::max(-1, boundsMaxFail(g_kbounds, p1_key));
        int side2_fail_until = std::max(-1, boundsMaxFail(g_kbounds, p2_key));
        int side1_true_from = std::max(-1, boundsMinSuccess(g_kbounds, p1_key));
        int side2_true_from = std::max(-1, boundsMinSuccess(g_kbounds, p2_key));
        if (side1_true_from < 0) side1_true_from = INT_MAX;
        if (side2_true_from < 0) side2_true_from = INT_MAX;

        for (int k1 = min_k1; k1 <= max_k1; ) {
            if (k1 <= side1_fail_until) {
                k1 = side1_fail_until + 1;
                continue;
            }

            int k2 = k - k1;
            if (k2 <= side2_fail_until) {
                return false;
            }

            const bool s1_guaranteed = k1 >= side1_true_from;
            const bool s2_guaranteed = k2 >= side2_true_from;
            if (s1_guaranteed && s2_guaranteed) {
                return true;
            }

            if (s1_guaranteed) {
                bool right_known = false;
                bool right_ok = false;
                right_known = tryBoundsPrune(g_kbounds, p2_key, k2, right_ok);
                if (!right_known) {
                    right_ok = FlipDistTree(A2, B2, k2);
                }
                if (right_ok) {
                    return true;
                }
                side2_fail_until = std::max(side2_fail_until,
                                            boundsMaxFail(g_kbounds, p2_key));
                ++k1;
                continue;
            }

            if (s2_guaranteed) {
                bool left_known = false;
                bool left_ok = false;
                left_known = tryBoundsPrune(g_kbounds, p1_key, k1, left_ok);
                if (!left_known) {
                    left_ok = FlipDistTree(A1, B1, k1);
                }
                if (left_ok) {
                    return true;
                }
                side1_fail_until = std::max(side1_fail_until,
                                            boundsMaxFail(g_kbounds, p1_key));
                ++k1;
                continue;
            }

            if (lb1 <= lb2) {
                bool left_known = false;
                bool left_ok = false;
                left_known = tryBoundsPrune(g_kbounds, p1_key, k1, left_ok);
                if (!left_known) {
                    left_ok = FlipDistTree(A1, B1, k1);
                }
                if (!left_ok) {
                    side1_fail_until = std::max(side1_fail_until,
                                                boundsMaxFail(g_kbounds, p1_key));
                    k1 = std::max(k1, side1_fail_until + 1);
                    continue;
                }
                side1_true_from = std::min(side1_true_from, k1);

                bool right_known = false;
                bool right_ok = false;
                right_known = tryBoundsPrune(g_kbounds, p2_key, k2, right_ok);
                if (!right_known) {
                    right_ok = FlipDistTree(A2, B2, k2);
                }
                if (right_ok) {
                    return true;
                }
                side2_fail_until = std::max(side2_fail_until,
                                            boundsMaxFail(g_kbounds, p2_key));
                ++k1;
            } else {
                bool right_known = false;
                bool right_ok = false;
                right_known = tryBoundsPrune(g_kbounds, p2_key, k2, right_ok);
                if (!right_known) {
                    right_ok = FlipDistTree(A2, B2, k2);
                }
                if (!right_ok) {
                    side2_fail_until = std::max(side2_fail_until,
                                                boundsMaxFail(g_kbounds, p2_key));
                    ++k1;
                    continue;
                }
                side2_true_from = std::min(side2_true_from, k2);

                bool left_known = false;
                bool left_ok = false;
                left_known = tryBoundsPrune(g_kbounds, p1_key, k1, left_ok);
                if (!left_known) {
                    left_ok = FlipDistTree(A1, B1, k1);
                }
                if (left_ok) {
                    return true;
                }
                side1_fail_until = std::max(side1_fail_until,
                                            boundsMaxFail(g_kbounds, p1_key));
                ++k1;
            }
        }
    }

    return false;
}

// Definition: Enumerate all independent edge subsets
// Parameters: edges: full edge list; index: current position; current: working subset; result: output list
// Returns: nothing; result is appended in-place
// Errors: none; may be expensive for large edge sets
void generateAllIndependentSubsets(const std::vector<std::pair<int,int>>& edges, int index,
                                   std::vector<std::pair<int,int>>& current,
                                   std::vector<std::vector<std::pair<int,int>>>& result) {
    if (index == edges.size()) {
        result.push_back(current);
        return;
    }

    // Choice 1: exclude current edge
    generateAllIndependentSubsets(edges, index + 1, current, result);

    // Choice 2: include current edge if it doesn't conflict
    bool canInclude = true;
    for (const auto& e : current) {
        if (areAdjacent(e, edges[index])) {
            canInclude = false;
            break;
        }
    }

    if (canInclude) {
        current.push_back(edges[index]);
        generateAllIndependentSubsets(edges, index + 1, current, result);
        current.pop_back();
    }
}

// Definition: Collect all edges incident to a node (parent + children)
// Parameters: T: tree; node: node id
// Returns: list of incident (parent, child) edges
// Errors: returns empty on invalid nodes or tree access failures
std::vector<std::pair<int,int>> getIncidentEdges(const VectorRangeTreeMap& T, int node) {
    std::vector<std::pair<int,int>> incident;

    try {
        if (!T.isOriginal(node)) return incident;

        int left = T.getLeftChild(node);
        int right = T.getRightChild(node);
        if (left >= 0 && T.isOriginal(left)) {
            incident.emplace_back(node, left);
        }
        if (right >= 0 && T.isOriginal(right)) {
            incident.emplace_back(node, right);
        }

        int parent = T.getParent(node);
        if (parent >= 0 && T.isOriginal(parent)) {
            incident.emplace_back(parent, node);
        }
    } catch (...) {
        // Return empty on error
    }

    return incident;
}

// Definition: Split S into pairs that fall entirely in T1 vs T2 node sets
// Parameters: S: edge-pair list; T1/T2: subtree partitions
// Returns: {S1, S2} containing pairs fully inside each subtree
// Errors: cross-partition pairs are dropped
std::pair<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>,
        std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>>
partitionS(const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S,
           const VectorRangeTreeMap& T1, const VectorRangeTreeMap& T2) {

    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S1, S2;
    if (S.empty()) {
        return {S1, S2};
    }
    if (T1.original_nodes.empty() || T2.original_nodes.empty()) {
        return {S1, S2};
    }

    S1.reserve(S.size());
    S2.reserve(S.size());

    const int max_node = std::max(T1.max_node_value, T2.max_node_value);
    if (max_node <= 0) {
        return {S1, S2};
    }

    std::vector<std::uint8_t> in_nodes_1(static_cast<std::size_t>(max_node + 1), 0);
    std::vector<std::uint8_t> in_nodes_2(static_cast<std::size_t>(max_node + 1), 0);
    for (int node : T1.original_nodes) {
        if (node >= 0 && node <= max_node) {
            in_nodes_1[static_cast<std::size_t>(node)] = 1;
        }
    }
    for (int node : T2.original_nodes) {
        if (node >= 0 && node <= max_node) {
            in_nodes_2[static_cast<std::size_t>(node)] = 1;
        }
    }

    for (const auto& edgePair : S) {
        auto& edge1 = edgePair.first;
        auto& edge2 = edgePair.second;

        // Check if both edges of the pair belong to T1
        bool edge1_in_T1 = edge1.first >= 0 && edge1.first <= max_node &&
                           edge1.second >= 0 && edge1.second <= max_node &&
                           in_nodes_1[static_cast<std::size_t>(edge1.first)] &&
                           in_nodes_1[static_cast<std::size_t>(edge1.second)];
        bool edge2_in_T1 = edge2.first >= 0 && edge2.first <= max_node &&
                           edge2.second >= 0 && edge2.second <= max_node &&
                           in_nodes_1[static_cast<std::size_t>(edge2.first)] &&
                           in_nodes_1[static_cast<std::size_t>(edge2.second)];

        // Check if both edges of the pair belong to T2
        bool edge1_in_T2 = edge1.first >= 0 && edge1.first <= max_node &&
                           edge1.second >= 0 && edge1.second <= max_node &&
                           in_nodes_2[static_cast<std::size_t>(edge1.first)] &&
                           in_nodes_2[static_cast<std::size_t>(edge1.second)];
        bool edge2_in_T2 = edge2.first >= 0 && edge2.first <= max_node &&
                           edge2.second >= 0 && edge2.second <= max_node &&
                           in_nodes_2[static_cast<std::size_t>(edge2.first)] &&
                           in_nodes_2[static_cast<std::size_t>(edge2.second)];

        if (edge1_in_T1 && edge2_in_T1) {
            S1.push_back(edgePair);
        } else if (edge1_in_T2 && edge2_in_T2) {
            S2.push_back(edgePair);
        }
        // If edge pair spans both partitions, we could assign to both or neither
        // For simplicity, we'll ignore cross-partition pairs
    }

    return {S1, S2};
}

bool partitionNodeSetsMatchByChildRange(const VectorRangeTreeMap& left,
                                       const VectorRangeTreeMap& right,
                                       const std::pair<int, int>& child_range) {
    try {
        int child = -1;
        for (int v : left.original_nodes) {
            if (!left.isOriginal(v)) continue;
            auto r = left.getRange(v);
            if (r.first == child_range.first && r.second == child_range.second) {
                child = v;
                break;
            }
        }
        if (child < 0) return false;
        int parent = left.getParent(child);
        if (parent < 0 || !left.isOriginal(parent)) {
            return false;
        }

        if (left.original_nodes != right.original_nodes) {
            return false;
        }
        for (int node : left.original_inorder) {
            if (!left.isOriginal(node) || !right.isOriginal(node)) {
                return false;
            }
            int pos = -1;
            if (node >= 0 && node < static_cast<int>(left.position_in_inorder_fast.size())) {
                pos = left.position_in_inorder_fast[node];
            } else {
                auto it = left.position_in_inorder.find(node);
                if (it != left.position_in_inorder.end()) {
                    pos = it->second;
                }
            }
            const bool in_left_child = pos >= child_range.first && pos < child_range.second;
            auto right_range = right.getRange(node);
            const bool in_right_child = right_range.first >= child_range.first &&
                                        right_range.second <= child_range.second;
            if (in_left_child != in_right_child) {
                return false;
            }
        }
        return true;
    } catch (...) {
        return false;
    }
}

void appendPartnerPairsFromSingleDiagonalProfiled(
    const VectorRangeTreeMap& T,
    const std::pair<int, int>& diagonal,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& out_pairs,
    PrefixPairgenRequestStats* stats) {
    PrefixPairgenRequestStats local;
    if (!stats) {
        stats = &local;
    }

    std::unordered_map<int, std::vector<EndpointEntry>> endpointMap;
    buildEndpointIndex(T, endpointMap);
    int vcount = vertexCount(T);
    stats->vcount = vcount;

    int L = diagonal.first;
    int R = diagonal.second;

    auto p1 = pickNeighborPair(L, R, vcount, endpointMap);
    if (!(p1.first.first < 0 && p1.second.first < 0)) {
        out_pairs.push_back(p1);
    }

    auto p2 = pickNeighborPair(R, L, vcount, endpointMap);
    if (!(p2.first.first < 0 && p2.second.first < 0)) {
        out_pairs.push_back(p2);
    }

    auto itL = endpointMap.find(L);
    auto itR = endpointMap.find(R);
    stats->bucket_size_l = (itL == endpointMap.end()) ? 0 : static_cast<int>(itL->second.size());
    stats->bucket_size_r = (itR == endpointMap.end()) ? 0 : static_cast<int>(itR->second.size());
    stats->emitted_pairs = static_cast<int>(out_pairs.size());
}

// Definition: Build independent edge subsets from S
// Parameters: S: list of edge pairs to branch on
// Returns: vector of independent edge sets
// Errors: none; may be exponential in size of S
std::vector<std::vector<std::pair<int,int>>> generateIndependentSubsetsFromS(
        const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {

    // First, collect all unique edges from S
    std::set<std::pair<int,int>> allEdges;
    for (const auto& edgePair : S) {
        allEdges.insert(edgePair.first);
        allEdges.insert(edgePair.second);
    }

    // Convert to vector for subset generation
    std::vector<std::pair<int,int>> edgeVector(allEdges.begin(), allEdges.end());

    // Generate independent subsets
    std::vector<std::vector<std::pair<int,int>>> independentSubsets;
    std::vector<std::pair<int,int>> current;

    generateAllIndependentSubsets(edgeVector, 0, current, independentSubsets);

    return independentSubsets;
}

void appendPartnerPairsFromDiagonals(
    const VectorRangeTreeMap &T,
    const std::vector<std::pair<int, int>> &diagonals,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &out_pairs) {

    std::unordered_map<int, std::vector<EndpointEntry>> endpointMap;
    buildEndpointIndex(T, endpointMap);
    int vcount = vertexCount(T);

    for (const auto &diag : diagonals) {
        int L = diag.first;
        int R = diag.second;

        auto p1 = pickNeighborPair(L, R, vcount, endpointMap);
        if (!(p1.first.first < 0 && p1.second.first < 0)) {
            out_pairs.emplace_back(p1);
        }

        auto p2 = pickNeighborPair(R, L, vcount, endpointMap);
        if (!(p2.first.first < 0 && p2.second.first < 0)) {
            out_pairs.emplace_back(p2);
        }
    }
}
