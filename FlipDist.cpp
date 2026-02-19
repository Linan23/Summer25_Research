#include "bf_bst.h"
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
#include <sstream>
#include <cstdlib>

// Forward declarations
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k);
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I);
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S);

// Debug flag
const bool DEBUG = false;  // Turn off debug for cleaner testing output

struct ProfileStats {
    bool enabled = false;
    long long calls_flipdist = 0;
    long long calls_tdi = 0;
    long long calls_tds = 0;
    long long calls_partition = 0;
    long long free_edge_hits = 0;
    long long free_edge_misses = 0;
    long long s_branch_calls = 0;
    long long s_empty_branch_calls = 0;
    long long independent_subsets_initial = 0;
    long long independent_subsets_s = 0;
    long long max_s_size = 0;
    long long max_i_size = 0;
    double time_flipdist_ms = 0.0;
    double time_tdi_ms = 0.0;
    double time_tds_ms = 0.0;
    double time_partition_ms = 0.0;
    std::chrono::steady_clock::time_point start_time;
    int abort_ms = -1;
};

static ProfileStats g_profile;

struct Key128 {
    uint64_t hi = 0;
    uint64_t lo = 0;
    bool operator==(const Key128& o) const {
        return hi == o.hi && lo == o.lo;
    }
};

struct Key128Hash {
    size_t operator()(const Key128& k) const {
        uint64_t h = k.hi;
        h ^= (k.lo + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        return static_cast<size_t>(h);
    }
};

struct TreeFingerprint {
    uint64_t h1 = 0;
    uint64_t h2 = 0;
};

static std::unordered_map<Key128, bool, Key128Hash> g_memo_flipdist;
static std::unordered_map<Key128, bool, Key128Hash> g_memo_tdi;
static std::unordered_map<Key128, bool, Key128Hash> g_memo_tds;
struct KBounds {
    int max_fail = -1;
    int min_success = -1;
};
static std::unordered_map<Key128, KBounds, Key128Hash> g_kbounds;
static std::unordered_map<Key128, KBounds, Key128Hash> g_tds_bounds;
static std::unordered_map<Key128, KBounds, Key128Hash> g_tdi_bounds;

static uint64_t fnvMix(uint64_t h, uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
    return h;
}

static uint64_t edgeKey(const std::pair<int,int>& e) {
    int a = std::min(e.first, e.second);
    int b = std::max(e.first, e.second);
    return (static_cast<uint64_t>(a) << 32) ^ static_cast<uint32_t>(b);
}

static uint64_t hashEdgeSet(const std::vector<std::pair<int,int>>& edges) {
    std::vector<uint64_t> keys;
    keys.reserve(edges.size());
    for (const auto& e : edges) keys.push_back(edgeKey(e));
    std::sort(keys.begin(), keys.end());
    uint64_t h = 1469598103934665603ULL;
    for (uint64_t v : keys) {
        h = fnvMix(h, v);
    }
    return h;
}

static uint64_t hashPairSet(const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& pairs) {
    std::vector<std::pair<uint64_t,uint64_t>> keys;
    keys.reserve(pairs.size());
    for (const auto& p : pairs) {
        uint64_t k1 = edgeKey(p.first);
        uint64_t k2 = edgeKey(p.second);
        if (k2 < k1) std::swap(k1, k2);
        keys.emplace_back(k1, k2);
    }
    std::sort(keys.begin(), keys.end());
    uint64_t h = 1469598103934665603ULL;
    for (const auto& kv : keys) {
        h = fnvMix(h, kv.first);
        h = fnvMix(h, kv.second);
    }
    return h;
}

static void dedupeSPairs(std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& pairs) {
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> out;
    out.reserve(pairs.size());
    std::unordered_set<Key128, Key128Hash> seen;
    for (const auto& p : pairs) {
        uint64_t k1 = edgeKey(p.first);
        uint64_t k2 = edgeKey(p.second);
        if (k2 < k1) std::swap(k1, k2);
        Key128 key{k1, k2};
        if (seen.insert(key).second) {
            out.push_back(p);
        }
    }
    pairs.swap(out);
}

static void mixKey(Key128& k, uint64_t v) {
    k.hi = fnvMix(k.hi, v);
    k.lo = fnvMix(k.lo, v ^ 0x9e3779b97f4a7c15ULL);
}

static void mixFingerprint(TreeFingerprint& fp, uint64_t v) {
    fp.h1 = fnvMix(fp.h1, v);
    fp.h2 = fnvMix(fp.h2, v ^ 0x9e3779b97f4a7c15ULL);
}

static TreeFingerprint treeFingerprint(const VectorRangeTreeMap& T) {
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());
    TreeFingerprint out{1469598103934665603ULL, 1099511628211ULL};
    mixFingerprint(out, static_cast<uint64_t>(nodes.size()));
    mixFingerprint(out, static_cast<uint64_t>(T.root + 1));
    for (int v : nodes) {
        auto r = T.getRange(v);
        int L = T.getLeftChild(v);
        int R = T.getRightChild(v);
        int P = T.getParent(v);
        mixFingerprint(out, static_cast<uint64_t>(v + 1));
        mixFingerprint(out, (static_cast<uint64_t>(static_cast<uint32_t>(r.first)) << 32) ^
                            static_cast<uint32_t>(r.second));
        mixFingerprint(out, (static_cast<uint64_t>(static_cast<uint32_t>(L + 2)) << 32) ^
                            static_cast<uint32_t>(R + 2));
        mixFingerprint(out, static_cast<uint64_t>(P + 2));
    }
    return out;
}

static Key128 makeKeyPair(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 1);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    return key;
}

static Key128 makeKeyFlip(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 2);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, static_cast<uint64_t>(k));
    return key;
}

static Key128 makeKeyI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
                       const std::vector<std::pair<int,int>>& I) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    uint64_t h = hashEdgeSet(I);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 3);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(I.size()));
    mixKey(key, static_cast<uint64_t>(k));
    return key;
}

static Key128 makeKeyS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
                       const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    uint64_t h = hashPairSet(S);
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_end);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 4);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(S.size()));
    mixKey(key, static_cast<uint64_t>(k));
    return key;
}

static Key128 makeKeySBase(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end,
                           const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    uint64_t h = hashPairSet(S);
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_end);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 5);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(S.size()));
    return key;
}

static Key128 makeKeyIBase(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final,
                           const std::vector<std::pair<int,int>>& I) {
    uint64_t h = hashEdgeSet(I);
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 6);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(I.size()));
    return key;
}

static bool tryBoundsPrune(const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                           const Key128& key, int k, bool& value_out) {
    auto it = bounds_map.find(key);
    if (it == bounds_map.end()) return false;
    const auto& b = it->second;
    if (b.min_success >= 0 && k >= b.min_success) {
        value_out = true;
        return true;
    }
    if (k <= b.max_fail) {
        value_out = false;
        return true;
    }
    return false;
}

static void updateBounds(std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                         const Key128& key, int k, bool value) {
    auto& b = bounds_map[key];
    if (value) {
        if (b.min_success < 0 || k < b.min_success) b.min_success = k;
        if (b.max_fail >= b.min_success) b.max_fail = b.min_success - 1;
    } else {
        if (k > b.max_fail) b.max_fail = k;
    }
}

static int requiredBudgetFromBounds(const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                    const Key128& key) {
    auto it = bounds_map.find(key);
    if (it == bounds_map.end()) return 0;
    return std::max(0, it->second.max_fail + 1);
}

static void resetMemo() {
    g_memo_flipdist.clear();
    g_memo_tdi.clear();
    g_memo_tds.clear();
    g_kbounds.clear();
    g_tds_bounds.clear();
    g_tdi_bounds.clear();
}

static void initProfile() {
    if (g_profile.enabled) return;
    const char *env = std::getenv("FLIPDIST_PROFILE");
    if (env && std::string(env) == "1") {
        g_profile.enabled = true;
        g_profile.start_time = std::chrono::steady_clock::now();
        const char *abort_env = std::getenv("FLIPDIST_PROFILE_ABORT_MS");
        if (abort_env) {
            g_profile.abort_ms = std::atoi(abort_env);
        }
    }
}

static bool profileAbortRequested() {
    if (!g_profile.enabled || g_profile.abort_ms <= 0) return false;
    auto now = std::chrono::steady_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(now - g_profile.start_time).count();
    return elapsed_ms > g_profile.abort_ms;
}

struct ScopedTimer {
    double *acc = nullptr;
    std::chrono::steady_clock::time_point t0;
    ScopedTimer(double *a) : acc(a) {
        if (g_profile.enabled && acc) {
            t0 = std::chrono::steady_clock::now();
        }
    }
    ~ScopedTimer() {
        if (g_profile.enabled && acc) {
            auto t1 = std::chrono::steady_clock::now();
            *acc += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
    }
};

// Definition: Emit a debug line when DEBUG is enabled
// Parameters: msg: message to print
// Returns: nothing
// Errors: none
void debugPrint(const std::string& msg) {
    if (DEBUG) {
        std::cout << "[DEBUG] " << msg << std::endl;
    }
}

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
    return getInternalEdges(T).size();
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

// Definition: Build a deep copy via preorder/inorder reconstruction
// Parameters: T: tree to copy
// Returns: copied tree, or an empty tree on failure
// Errors: returns empty on invalid trees or reconstruction failure
VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap& T) {
    VectorRangeTreeMap copy;
    try {
        if (T.root < 0 || !T.isOriginal(T.root) || T.original_nodes.empty()) {
            return copy;
        }

        std::vector<int> preorder, inorder;

        std::function<void(int, std::vector<int>&)> buildPreorder = [&](int node, std::vector<int>& pre) {
            if (node < 0 || !T.isOriginal(node)) return;
            pre.push_back(node);
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left)) buildPreorder(left, pre);
            if (right >= 0 && T.isOriginal(right)) buildPreorder(right, pre);
        };

        std::function<void(int, std::vector<int>&)> buildInorder = [&](int node, std::vector<int>& in) {
            if (node < 0 || !T.isOriginal(node)) return;
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left)) buildInorder(left, in);
            in.push_back(node);
            if (right >= 0 && T.isOriginal(right)) buildInorder(right, in);
        };

        buildPreorder(T.root, preorder);
        buildInorder(T.root, inorder);

        if (preorder.size() == inorder.size() && !preorder.empty()) {
            copy.build(preorder, inorder);
        }
    } catch (...) {
        VectorRangeTreeMap empty;
        return empty;
    }
    return copy;
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
static bool tryCommonEdgeDecomposition(const VectorRangeTreeMap& T_init,
                                       const VectorRangeTreeMap& T_final,
                                       int k,
                                       bool &handled) {
    handled = false;
    if (T_init.original_nodes.empty() || T_final.original_nodes.empty()) {
        return false;
    }

    auto targetEdges = buildTargetSet(T_final);

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

    for (const auto &edge : commonEdges) {
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
        int lb1 = std::max(conf1, requiredBudgetFromBounds(g_kbounds, p1_key));
        int lb2 = std::max(conf2, requiredBudgetFromBounds(g_kbounds, p2_key));

        int min_k1 = std::max(0, lb1);
        int max_k1 = std::min(k, k - lb2);
        if (min_k1 > max_k1) {
            continue;
        }

        for (int k1 = min_k1; k1 <= max_k1; k1++) {
            int k2 = k - k1;

            bool left_known = false;
            bool left_ok = false;
            left_known = tryBoundsPrune(g_kbounds, p1_key, k1, left_ok);
            if (!left_known) {
                left_ok = FlipDistTree(A1, B1, k1);
            }
            if (!left_ok) continue;

            bool right_known = false;
            bool right_ok = false;
            right_known = tryBoundsPrune(g_kbounds, p2_key, k2, right_ok);
            if (!right_known) {
                right_ok = FlipDistTree(A2, B2, k2);
            }
            if (right_ok) return true;
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

    // Get node sets for each partition
    std::set<int> nodes1, nodes2;
    for (int node : T1.original_nodes) nodes1.insert(node);
    for (int node : T2.original_nodes) nodes2.insert(node);

    for (const auto& edgePair : S) {
        auto& edge1 = edgePair.first;
        auto& edge2 = edgePair.second;

        // Check if both edges of the pair belong to T1
        bool edge1_in_T1 = nodes1.count(edge1.first) && nodes1.count(edge1.second);
        bool edge2_in_T1 = nodes1.count(edge2.first) && nodes1.count(edge2.second);

        // Check if both edges of the pair belong to T2
        bool edge1_in_T2 = nodes2.count(edge1.first) && nodes2.count(edge1.second);
        bool edge2_in_T2 = nodes2.count(edge2.first) && nodes2.count(edge2.second);

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

// Definition: Decide if rotation distance <= k via Li–Xia branching
// Parameters: T_init: source tree; T_final: target tree; k: budget
// Returns: true if a sequence of <= k rotations exists, else false
// Errors: returns false on invalid inputs or internal copy/rotation failures
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k) {
    initProfile();
    if (g_profile.enabled) {
        g_profile.calls_flipdist++;
    }
    ScopedTimer _t(g_profile.enabled ? &g_profile.time_flipdist_ms : nullptr);
    if (profileAbortRequested()) {
        return false;
    }
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    const Key128 memo_key = makeKeyFlip(T_init, T_final, k);
    const Key128 pair_key = makeKeyPair(T_init, T_final);
    auto memo_it = g_memo_flipdist.find(memo_key);
    if (memo_it != g_memo_flipdist.end()) {
        return memo_it->second;
    }
    bool pair_bounds_value = false;
    if (tryBoundsPrune(g_kbounds, pair_key, k, pair_bounds_value)) {
        g_memo_flipdist[memo_key] = pair_bounds_value;
        return pair_bounds_value;
    }
    auto memoReturn = [&](bool value) {
        g_memo_flipdist[memo_key] = value;
        updateBounds(g_kbounds, pair_key, k, value);
        return value;
    };

    // BASE CASE: Check if trees are already identical
    if (TreesEqual(T_init, T_final)) {
        debugPrint("Trees already equal, returning true");
        return memoReturn(true);
    }

    bool handled_common = false;
    if (tryCommonEdgeDecomposition(T_init, T_final, k, handled_common)) {
        debugPrint("Common edge decomposition succeeded");
        return memoReturn(true);
    }
    if (handled_common) {
        debugPrint("Common edge decomposition failed");
        return memoReturn(false);
    }

    auto [hasFree, freeEdge] = findFreeEdge(T_init, T_final);
    if (hasFree) {
        debugPrint("FlipDistTree: free edge found, delegating to TreeDistS");
        return memoReturn(TreeDistS(T_init, T_final, k, {}));
    }

    debugPrint("FlipDistTree: budget k=" + std::to_string(k));

    auto edges = getInternalEdges(T_init);
    debugPrint("Found " + std::to_string(edges.size()) + " internal edges");

    std::vector<std::vector<std::pair<int,int>>> independentSubsets;
    std::vector<std::pair<int,int>> current;

    //                       branch on two choices: (1) include e in I; (2) exclude e from I"
    generateAllIndependentSubsets(edges, 0, current, independentSubsets);
    if (g_profile.enabled) {
        g_profile.independent_subsets_initial += (long long)independentSubsets.size();
    }
    debugPrint("Generated " + std::to_string(independentSubsets.size()) + " independent subsets (all)");

    //                       If FlipDistTree-I(T_init, T_final, |I|) returns True, then True"
    for (size_t i = 0; i < independentSubsets.size(); i++) {
        const auto& subset = independentSubsets[i];
        if (subset.empty()) continue;  // Skip empty subsets (I ≠ ∅ requirement)
        if (g_profile.enabled) {
            if ((long long)subset.size() > g_profile.max_i_size) {
                g_profile.max_i_size = (long long)subset.size();
            }
        }

        debugPrint("Trying subset " + std::to_string(i) + " with " + std::to_string(subset.size()) + " edges");

        try {
            // Call TreeDistI (which is FlipDistTree-I from pseudocode)
            if (TreeDistI(T_init, T_final, k, subset)) {
                debugPrint("Found solution with subset " + std::to_string(i));
                return memoReturn(true);
            }
        } catch (...) {
            debugPrint("Exception in TreeDistI for subset " + std::to_string(i));
            continue;
        }
    }

    debugPrint("No solution found, returning false");
    return memoReturn(false);
}

// Definition: Apply an independent set of rotations, then branch on S
// Parameters: T_init/T_final: trees; k: total budget; I: independent edges
// Returns: true if solvable within k, else false
// Errors: returns false on invalid edges or copy/rotation failures
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I) {
    initProfile();
    if (g_profile.enabled) {
        g_profile.calls_tdi++;
        if ((long long)I.size() > g_profile.max_i_size) {
            g_profile.max_i_size = (long long)I.size();
        }
    }
    ScopedTimer _t(g_profile.enabled ? &g_profile.time_tdi_ms : nullptr);
    if (profileAbortRequested()) {
        return false;
    }
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    const Key128 memo_key = makeKeyI(T_init, T_final, k, I);
    auto memo_it = g_memo_tdi.find(memo_key);
    if (memo_it != g_memo_tdi.end()) {
        return memo_it->second;
    }
    const Key128 bounds_key = makeKeyIBase(T_init, T_final, I);
    bool bounds_value = false;
    if (tryBoundsPrune(g_tdi_bounds, bounds_key, k, bounds_value)) {
        g_memo_tdi[memo_key] = bounds_value;
        return bounds_value;
    }
    auto memoReturn = [&](bool value) {
        g_memo_tdi[memo_key] = value;
        updateBounds(g_tdi_bounds, bounds_key, k, value);
        return value;
    };

    // Exact lower bound: each conflicting edge must be removed by at least one rotation
    int conflict_lb = countConflictEdges(T_init, T_final);
    if (conflict_lb > k) {
        return memoReturn(false);
    }

    int remaining_budget = k - (int)I.size();  // k - |I| from pseudocode

    if (remaining_budget < 0) {  // This covers φ(T_init) > k − |I| case
        debugPrint("TreeDistI: Not enough budget for |I| rotations");
        return memoReturn(false);
    }

    // Special handling when budget exactly equals |I| - try direct solution
    if (remaining_budget == 0) {
        VectorRangeTreeMap T_test = safeCopyTree(T_init);

        // Apply all rotations in I
        for (const auto& edge : I) {
            auto oriented = orientEdge(T_test, edge);
            int parent = oriented.first;
            int child = oriented.second;
            if (parent < 0 || child < 0) continue;

            try {
                if (T_test.getLeftChild(parent) == child) {
                    T_test.rotateRight(parent);
                } else if (T_test.getRightChild(parent) == child) {
                    T_test.rotateLeft(parent);
                }
            } catch (...) {
                continue;
            }
        }

        if (TreesEqual(T_test, T_final)) {
            debugPrint("TreeDistI: Solved with exactly |I| rotations");
            return memoReturn(true);
        } else {
            debugPrint("TreeDistI: No budget left and not solved");
            return memoReturn(false);
        }
    }

    debugPrint("TreeDistI: Proceeding with remaining_budget=" + std::to_string(remaining_budget));

    if (TreesEqual(T_init, T_final)) {
        debugPrint("TreeDistI: Trees already equal");
        return memoReturn(true);
    }

    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S;
    std::vector<std::pair<int,int>> newDiagonals;

    VectorRangeTreeMap T_bar = safeCopyTree(T_init);  // T̄_init from pseudocode
    if (T_bar.original_nodes.empty()) {
        debugPrint("TreeDistI: Failed to copy tree");
        return memoReturn(false);
    }

    for (const auto& edge : I) {
        auto oriented = orientEdge(T_bar, edge);
        int parent = oriented.first;
        int child = oriented.second;

        debugPrint("TreeDistI: Processing edge (" + std::to_string(parent) + "," + std::to_string(child) + ")");

        if (parent < 0 || child < 0) {
            debugPrint("TreeDistI: Invalid edge, skipping");
            continue;
        }

        int u, v;  // The nodes u,v from pseudocode step 2.2

        if (T_bar.getLeftChild(parent) == child) {
            u = child;   // After rotation, child becomes parent
            v = parent;  // After rotation, parent becomes child
            T_bar.rotateRight(parent);
        } else if (T_bar.getRightChild(parent) == child) {
            u = child;   // After rotation, child becomes parent
            v = parent;  // After rotation, parent becomes child
            T_bar.rotateLeft(parent);
        } else {
            continue;
        }

        debugPrint("TreeDistI: Applied rotation, new edge ē connects " + std::to_string(u) + " and " + std::to_string(v));

        auto cr = T_bar.getRange(v);
        newDiagonals.emplace_back(cr.first, cr.second);
    }

    std::unordered_map<int, std::vector<EndpointEntry>> endpointMap;
    buildEndpointIndex(T_bar, endpointMap);
    int vcount = vertexCount(T_bar);

    for (const auto& diag : newDiagonals) {
        int L = diag.first;
        int R = diag.second;

        auto p1 = pickNeighborPair(L, R, vcount, endpointMap);
        if (!(p1.first.first < 0 && p1.second.first < 0)) {
            S.emplace_back(p1);
        }

        auto p2 = pickNeighborPair(R, L, vcount, endpointMap);
        if (!(p2.first.first < 0 && p2.second.first < 0)) {
            S.emplace_back(p2);
        }
    }

    int post_i_conflict_lb = countConflictEdges(T_bar, T_final);
    if (post_i_conflict_lb > k - (int)I.size()) {
        return memoReturn(false);
    }

    return memoReturn(TreeDistS(T_bar, T_final, k - (int)I.size(), S));
}

// Definition: Resolve free-edge/partition cases, else branch on S pairs
// Parameters: T_init/T_end: trees; k: budget; S: edge-pair set
// Returns: true if solvable within k, else false
// Errors: returns false on invalid inputs or internal failures
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    initProfile();
    if (g_profile.enabled) {
        g_profile.calls_tds++;
        if ((long long)S.size() > g_profile.max_s_size) {
            g_profile.max_s_size = (long long)S.size();
        }
    }
    ScopedTimer _t(g_profile.enabled ? &g_profile.time_tds_ms : nullptr);
    if (profileAbortRequested()) {
        return false;
    }
    debugPrint("Entering TreeDistS with k=" + std::to_string(k) + ", |S|=" + std::to_string(S.size()));

    const Key128 memo_key = makeKeyS(T_init, T_end, k, S);
    auto memo_it = g_memo_tds.find(memo_key);
    if (memo_it != g_memo_tds.end()) {
        return memo_it->second;
    }
    const Key128 bounds_key = makeKeySBase(T_init, T_end, S);
    bool bounds_value = false;
    if (tryBoundsPrune(g_tds_bounds, bounds_key, k, bounds_value)) {
        g_memo_tds[memo_key] = bounds_value;
        return bounds_value;
    }
    auto memoReturn = [&](bool value) {
        g_memo_tds[memo_key] = value;
        updateBounds(g_tds_bounds, bounds_key, k, value);
        return value;
    };

    // Base case: trees already equal
    if (TreesEqual(T_init, T_end)) {
        debugPrint("TreeDistS: Trees already equal");
        return memoReturn(true);
    }

    if (k < 0) {
        debugPrint("TreeDistS: Negative budget");
        return memoReturn(false);
    }
    if (k == 0) {
        debugPrint("TreeDistS: No budget and trees differ");
        return memoReturn(false);
    }

    int conflict_lb = countConflictEdges(T_init, T_end);
    if (conflict_lb > k) {
        return memoReturn(false);
    }

    VectorRangeTreeMap T_work = safeCopyTree(T_init);
    int k_work = k;
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_work = S;
    dedupeSPairs(S_work);
    struct UnfoldedStateMemoEntry {
        Key128 exact_key;
        Key128 bounds_key;
        int budget = 0;
    };
    std::vector<UnfoldedStateMemoEntry> unfolded_bounds_entries;

    auto commitResult = [&](bool value) {
        for (const auto& entry : unfolded_bounds_entries) {
            g_memo_tds[entry.exact_key] = value;
            updateBounds(g_tds_bounds, entry.bounds_key, entry.budget, value);
        }
        return memoReturn(value);
    };

    while (k_work > 0) {
        int work_conflict_lb = countConflictEdges(T_work, T_end);
        if (work_conflict_lb > k_work) {
            return commitResult(false);
        }

        const Key128 state_key = makeKeyS(T_work, T_end, k_work, S_work);
        auto state_it = g_memo_tds.find(state_key);
        if (state_it != g_memo_tds.end()) {
            return commitResult(state_it->second);
        }
        const Key128 state_bounds_key = makeKeySBase(T_work, T_end, S_work);
        bool state_bounds_value = false;
        if (tryBoundsPrune(g_tds_bounds, state_bounds_key, k_work, state_bounds_value)) {
            return commitResult(state_bounds_value);
        }
        unfolded_bounds_entries.push_back({state_key, state_bounds_key, k_work});

        auto [hasFree, freeEdge] = findFreeEdge(T_work, T_end);
        if (g_profile.enabled) {
            if (hasFree) g_profile.free_edge_hits++;
            else g_profile.free_edge_misses++;
        }
        if (!hasFree) {
            break;
        }

        debugPrint("TreeDistS: Found free edge (" + std::to_string(freeEdge.first) + "," + std::to_string(freeEdge.second) + ")");

        try {
            int parent = freeEdge.first;
            int child = freeEdge.second;

            std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_filtered;
            S_filtered.reserve(S_work.size() + 2);
            for (const auto& pair : S_work) {
                if (!(pair.first == freeEdge || pair.second == freeEdge)) {
                    S_filtered.push_back(pair);
                }
            }
            debugPrint("TreeDistS: Filtered S from " + std::to_string(S_work.size()) + " to " + std::to_string(S_filtered.size()) + " pairs");

            VectorRangeTreeMap T_bar = safeCopyTree(T_work);
            int u, v;

            if (T_bar.getLeftChild(parent) == child) {
                T_bar.rotateRight(parent);
                u = child;
                v = parent;
            } else if (T_bar.getRightChild(parent) == child) {
                T_bar.rotateLeft(parent);
                u = child;
                v = parent;
            } else {
                return commitResult(false);
            }

            if (TreesEqual(T_bar, T_end)) {
                debugPrint("TreeDistS: Solved with free edge rotation");
                return commitResult(true);
            }

            auto cr = T_bar.getRange(v);
            int L = cr.first;
            int R = cr.second;

            std::unordered_map<int, std::vector<EndpointEntry>> endpointMap;
            buildEndpointIndex(T_bar, endpointMap);
            int vcount = vertexCount(T_bar);

            auto p1 = pickNeighborPair(L, R, vcount, endpointMap);
            if (!(p1.first.first < 0 && p1.second.first < 0)) {
                S_filtered.emplace_back(p1);
            }
            auto p2 = pickNeighborPair(R, L, vcount, endpointMap);
            if (!(p2.first.first < 0 && p2.second.first < 0)) {
                S_filtered.emplace_back(p2);
            }
            dedupeSPairs(S_filtered);

            auto parent_range = T_bar.getRange(u);
            auto child_range = T_bar.getRange(v);

            try {
                if (g_profile.enabled) {
                    g_profile.calls_partition++;
                }
                ScopedTimer _tp(g_profile.enabled ? &g_profile.time_partition_ms : nullptr);
                auto [T_bar1, T_bar2] = VectorRangeTreeMap::partitionAlongEdge(T_bar, parent_range, child_range);
                auto [T_end1, T_end2] = VectorRangeTreeMap::partitionAlongEdge(T_end, parent_range, child_range);

                if (T_bar1.original_nodes != T_end1.original_nodes ||
                    T_bar2.original_nodes != T_end2.original_nodes) {
                    debugPrint("TreeDistS: Partition mismatch, continuing free-edge contraction");
                    T_work = std::move(T_bar);
                    S_work = std::move(S_filtered);
                    --k_work;
                    continue;
                }

                auto [S1, S2] = partitionS(S_filtered, T_bar1, T_bar2);
                int n1 = countInternalEdges(T_bar1);
                int n2 = countInternalEdges(T_bar2);
                int conf1 = countConflictEdges(T_bar1, T_end1);
                int conf2 = countConflictEdges(T_bar2, T_end2);

                if (n1 == 0) {
                    if (conf2 > k_work - 1) return commitResult(false);
                    return commitResult(TreeDistS(T_bar2, T_end2, k_work - 1, S2));
                }

                if (n2 == 0) {
                    if (conf1 > k_work - 1) return commitResult(false);
                    return commitResult(TreeDistS(T_bar1, T_end1, k_work - 1, S1));
                }

                const Key128 s1_bounds_key = makeKeySBase(T_bar1, T_end1, S1);
                const Key128 s2_bounds_key = makeKeySBase(T_bar2, T_end2, S2);
                const int lb1 = std::max(requiredBudgetFromBounds(g_tds_bounds, s1_bounds_key), conf1);
                const int lb2 = std::max(requiredBudgetFromBounds(g_tds_bounds, s2_bounds_key), conf2);

                int min_k1 = std::max(0, lb1);
                int max_k1 = std::min(k_work - 1, k_work - 1 - lb2);
                if (min_k1 > max_k1) {
                    return commitResult(false);
                }

                for (int k1 = min_k1; k1 <= max_k1; k1++) {
                    int k2 = k_work - 1 - k1;

                    bool s1_known = false;
                    bool s1_ok = false;
                    s1_known = tryBoundsPrune(g_tds_bounds, s1_bounds_key, k1, s1_ok);
                    if (!s1_known) {
                        s1_ok = TreeDistS(T_bar1, T_end1, k1, S1);
                    }
                    if (!s1_ok) {
                        continue;
                    }

                    bool s2_known = false;
                    bool s2_ok = false;
                    s2_known = tryBoundsPrune(g_tds_bounds, s2_bounds_key, k2, s2_ok);
                    if (!s2_known) {
                        s2_ok = TreeDistS(T_bar2, T_end2, k2, S2);
                    }
                    if (s2_ok) {
                        return commitResult(true);
                    }
                }
                return commitResult(false);
            } catch (...) {
                debugPrint("TreeDistS: Partitioning failed, continuing free-edge contraction");
                T_work = std::move(T_bar);
                S_work = std::move(S_filtered);
                --k_work;
                continue;
            }
        } catch (...) {
            debugPrint("TreeDistS: Exception during free edge handling");
            return commitResult(false);
        }
    }

    //                     For each nonempty independent subset I ⊆ ⋃ S (no two edges in I share a node):"
    debugPrint("TreeDistS: No free edge, implementing S branching (step 2)");

    if (k_work <= 0) {
        debugPrint("TreeDistS: No budget left");
        return commitResult(false);
    }

    if (S_work.empty()) {
        // No constraints from S, try any rotation
        if (g_profile.enabled) {
            g_profile.s_empty_branch_calls++;
        }
        debugPrint("TreeDistS: S is empty, trying any rotation");
        const int current_conflicts = countConflictEdges(T_work, T_end);
        const bool must_drop_conflict_now = (k_work == current_conflicts);
        const int child_k = k_work - 1;
        const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> empty_S;
        std::unordered_set<Key128, Key128Hash> seen_child_states;
        auto edges = getInternalEdges(T_work);
        struct RotationCandidate {
            std::pair<int,int> edge;
            int next_conflicts = INT_MAX;
            int free_hint = 0;
            Key128 child_exact_key{};
            Key128 child_bounds_key{};
        };
        std::vector<RotationCandidate> candidates;
        candidates.reserve(edges.size());
        std::vector<RotationCandidate> plateau_candidates;
        plateau_candidates.reserve(edges.size());
        auto byPromise = [](const RotationCandidate& a, const RotationCandidate& b) {
            if (a.next_conflicts != b.next_conflicts) {
                return a.next_conflicts < b.next_conflicts;
            }
            if (a.free_hint != b.free_hint) {
                return a.free_hint > b.free_hint;
            }
            if (a.edge.first != b.edge.first) {
                return a.edge.first < b.edge.first;
            }
            return a.edge.second < b.edge.second;
        };

        for (const auto& edge : edges) {
            try {
                VectorRangeTreeMap probe = safeCopyTree(T_work);
                int parent = edge.first;
                int child = edge.second;
                if (probe.getLeftChild(parent) == child) {
                    probe.rotateRight(parent);
                } else if (probe.getRightChild(parent) == child) {
                    probe.rotateLeft(parent);
                } else {
                    continue;
                }
                int next_conflicts = countConflictEdges(probe, T_end);
                if (next_conflicts > k_work - 1) {
                    continue;
                }
                Key128 child_exact_key = makeKeyS(probe, T_end, child_k, empty_S);
                if (!seen_child_states.insert(child_exact_key).second) {
                    continue;
                }
                auto child_memo_it = g_memo_tds.find(child_exact_key);
                if (child_memo_it != g_memo_tds.end()) {
                    if (child_memo_it->second) {
                        return commitResult(true);
                    }
                    continue;
                }
                Key128 child_bounds_key = makeKeySBase(probe, T_end, empty_S);
                bool known_tds = false;
                if (tryBoundsPrune(g_tds_bounds, child_bounds_key, child_k, known_tds)) {
                    if (known_tds) {
                        return commitResult(true);
                    }
                    continue;
                }
                // Exact lower-bound gate for this child state before recursing
                Key128 child_pair_key = makeKeyPair(probe, T_end);
                int child_lb = std::max(next_conflicts,
                                        std::max(requiredBudgetFromBounds(g_tds_bounds, child_bounds_key),
                                                 requiredBudgetFromBounds(g_kbounds, child_pair_key)));
                if (child_lb > child_k) {
                    continue;
                }
                auto [probe_has_free, _] = findFreeEdge(probe, T_end);
                RotationCandidate cand{edge, next_conflicts, probe_has_free ? 1 : 0,
                                       child_exact_key, child_bounds_key};
                if (must_drop_conflict_now && next_conflicts >= current_conflicts) {
                    continue;
                }
                if (next_conflicts < current_conflicts || cand.free_hint) {
                    candidates.push_back(cand);
                } else {
                    plateau_candidates.push_back(cand);
                }
            } catch (...) {
                continue;
            }
        }

        std::sort(candidates.begin(), candidates.end(), byPromise);
        std::sort(plateau_candidates.begin(), plateau_candidates.end(), byPromise);
        candidates.insert(candidates.end(), plateau_candidates.begin(), plateau_candidates.end());

        for (const auto& cand : candidates) {
            try {
                auto child_memo_it = g_memo_tds.find(cand.child_exact_key);
                if (child_memo_it != g_memo_tds.end()) {
                    if (child_memo_it->second) {
                        return commitResult(true);
                    }
                    continue;
                }
                bool known_tds = false;
                if (tryBoundsPrune(g_tds_bounds, cand.child_bounds_key, child_k, known_tds)) {
                    if (known_tds) {
                        return commitResult(true);
                    }
                    continue;
                }
                VectorRangeTreeMap T_test = safeCopyTree(T_work);
                int parent = cand.edge.first;
                int child = cand.edge.second;

                if (T_test.getLeftChild(parent) == child) {
                    T_test.rotateRight(parent);
                } else if (T_test.getRightChild(parent) == child) {
                    T_test.rotateLeft(parent);
                }

                if (TreeDistS(T_test, T_end, child_k, empty_S)) {
                    debugPrint("TreeDistS: Found solution with rotation");
                    return commitResult(true);
                }
            } catch (...) {
                continue;
            }
        }
    } else {
        //                       (a) include neither, (b) include eᵢ only, (c) include eᵢ′ only
        //                       (skip choices that pick a non-edge or conflict with independence)."
        if (g_profile.enabled) {
            g_profile.s_branch_calls++;
        }
        debugPrint("TreeDistS: Implementing complete S branching");

        std::vector<std::pair<int,int>> chosen;
        std::unordered_set<int> used_nodes;
        dedupeSPairs(S_work);

        bool use_mask_memo = (T_work.original_nodes.size() <= 63);
        std::unordered_map<int, int> node_to_bit;
        if (use_mask_memo) {
            int bit = 0;
            for (int node : T_work.original_nodes) {
                node_to_bit[node] = bit++;
            }
        }
        auto edgeMask = [&](const std::pair<int,int>& e, uint64_t& mask_out) -> bool {
            if (!use_mask_memo) return false;
            auto it1 = node_to_bit.find(e.first);
            auto it2 = node_to_bit.find(e.second);
            if (it1 == node_to_bit.end() || it2 == node_to_bit.end()) return false;
            if (it1->second >= 63 || it2->second >= 63) return false;
            mask_out = (1ULL << it1->second) | (1ULL << it2->second);
            return true;
        };
        uint64_t used_mask = 0;
        struct BranchStateKey {
            uint32_t idx = 0;
            uint64_t mask = 0;
            bool operator==(const BranchStateKey& o) const {
                return idx == o.idx && mask == o.mask;
            }
        };
        struct BranchStateHash {
            size_t operator()(const BranchStateKey& k) const {
                uint64_t h = (static_cast<uint64_t>(k.idx) << 32) ^ k.mask;
                h ^= (h >> 33);
                h *= 0xff51afd7ed558ccdULL;
                h ^= (h >> 33);
                return static_cast<size_t>(h);
            }
        };
        std::unordered_set<BranchStateKey, BranchStateHash> failed_states;

        auto canInclude = [&](const std::pair<int,int>& e) -> bool {
            auto oriented = orientEdge(T_work, e);
            if (oriented.first < 0 || oriented.second < 0) {
                return false;
            }
            if (use_mask_memo) {
                uint64_t m = 0;
                if (!edgeMask(e, m)) return false;
                return (used_mask & m) == 0;
            }
            return used_nodes.count(e.first) == 0 && used_nodes.count(e.second) == 0;
        };

        std::function<bool(size_t)> branchPairs = [&](size_t idx) -> bool {
            BranchStateKey state{static_cast<uint32_t>(idx), used_mask};
            if (use_mask_memo && failed_states.find(state) != failed_states.end()) {
                return false;
            }
            if ((int)chosen.size() > k_work) {
                if (use_mask_memo) failed_states.insert(state);
                return false;
            }
            if (idx == S_work.size()) {
                if (chosen.empty()) {
                    if (use_mask_memo) failed_states.insert(state);
                    return false;
                }
                const Key128 i_bounds_key = makeKeyIBase(T_work, T_end, chosen);
                bool i_bounds_value = false;
                if (tryBoundsPrune(g_tdi_bounds, i_bounds_key, k_work, i_bounds_value)) {
                    if (!i_bounds_value && use_mask_memo) failed_states.insert(state);
                    return i_bounds_value;
                }
                bool i_ok = TreeDistI(T_work, T_end, k_work, chosen);
                if (!i_ok && use_mask_memo) failed_states.insert(state);
                return i_ok;
            }

            // Option 1: include neither
            if (branchPairs(idx + 1)) return true;

            const auto& pair = S_work[idx];
            const auto& e1 = pair.first;
            const auto& e2 = pair.second;

            // Option 2: include e1
            if (canInclude(e1)) {
                uint64_t m1 = 0;
                bool has_m1 = edgeMask(e1, m1);
                chosen.push_back(e1);
                used_nodes.insert(e1.first);
                used_nodes.insert(e1.second);
                if (has_m1) used_mask |= m1;
                if (branchPairs(idx + 1)) return true;
                if (has_m1) used_mask &= ~m1;
                used_nodes.erase(e1.first);
                used_nodes.erase(e1.second);
                chosen.pop_back();
            }

            // Option 3: include e2 (avoid duplicate edge)
            if (!(e1 == e2) && canInclude(e2)) {
                uint64_t m2 = 0;
                bool has_m2 = edgeMask(e2, m2);
                chosen.push_back(e2);
                used_nodes.insert(e2.first);
                used_nodes.insert(e2.second);
                if (has_m2) used_mask |= m2;
                if (branchPairs(idx + 1)) return true;
                if (has_m2) used_mask &= ~m2;
                used_nodes.erase(e2.first);
                used_nodes.erase(e2.second);
                chosen.pop_back();
            }

            if (use_mask_memo) failed_states.insert(state);
            return false;
        };

        if (branchPairs(0)) {
            debugPrint("TreeDistS: Found solution in S branching");
            return commitResult(true);
        }
    }

    debugPrint("TreeDistS: No solution found");
    return commitResult(false);
}
// ============================================================================
// COMPREHENSIVE ACCURACY AND SCALABILITY TESTING SUITE
// ============================================================================


// Definition: Simple stopwatch for milliseconds and microseconds timing
// Parameters: none
// Returns: constructed timer
// Errors: none
class PerformanceTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
public:
    // Definition: Start or reset the timer
    // Parameters: none
    // Returns: nothing
    // Errors: none
    void start() { start_time = std::chrono::high_resolution_clock::now(); }

    // Definition: Elapsed time in microseconds since start
    // Parameters: none
    // Returns: elapsed microseconds
    // Errors: none
    long long getMicroseconds() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }

    // Definition: Elapsed time in milliseconds since start
    // Parameters: none
    // Returns: elapsed milliseconds
    // Errors: none
    double getMilliseconds() { return getMicroseconds() / 1000.0; }
};

// Definition: Deterministic tree shape generators for CLI inputs
// Parameters: none
// Returns: class with static generation helpers
// Errors: none
class TreeGenerator {
public:
    // Definition: Right-leaning chain tree
    // Parameters: n: number of nodes
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateRightChain(int n) {
        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++) {
            preorder.push_back(i);
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    // Definition: Left-leaning chain tree
    // Parameters: n: number of nodes
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateLeftChain(int n) {
        std::vector<int> preorder, inorder;
        for (int i = n; i >= 1; i--) {
            preorder.push_back(i);
        }
        for (int i = 1; i <= n; i++) {
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    // Definition: Balanced BST shape
    // Parameters: n: number of nodes
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateBalanced(int n) {
        if (n == 0) return {{}, {}};
        if (n == 1) return {{1}, {1}};
        if (n == 2) return {{2, 1}, {1, 2}};
        if (n == 3) return {{2, 1, 3}, {1, 2, 3}};
        if (n == 4) return {{3, 2, 1, 4}, {1, 2, 3, 4}};
        if (n == 5) return {{3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}};
        if (n == 6) return {{4, 2, 1, 3, 5, 6}, {1, 2, 3, 4, 5, 6}};
        if (n == 7) return {{4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7}};

        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++) inorder.push_back(i);

        std::function<void(int, int)> buildBalanced = [&](int start, int end) {
            if (start > end) return;
            int mid = (start + end) / 2;
            preorder.push_back(mid);
            buildBalanced(start, mid - 1);
            buildBalanced(mid + 1, end);
        };

        buildBalanced(1, n);
        return {preorder, inorder};
    }

    // Definition: Random BST shape from shuffled insertion order
    // Parameters: n: number of nodes; seed: shuffle seed
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateRandom(int n, int seed = 42) {
        if (n <= 0) return {{}, {}};

        struct Node {
            int v;
            Node *l = nullptr;
            Node *r = nullptr;
            explicit Node(int val) : v(val) {}
        };
        auto insert = [](Node *root, int val) -> Node * {
            if (!root) return new Node(val);
            Node *cur = root;
            while (true) {
                if (val < cur->v) {
                    if (cur->l) cur = cur->l;
                    else { cur->l = new Node(val); break; }
                } else {
                    if (cur->r) cur = cur->r;
                    else { cur->r = new Node(val); break; }
                }
            }
            return root;
        };
        auto preorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
            if (!node) return;
            out.push_back(node->v);
            self(node->l, out, self);
            self(node->r, out, self);
        };
        auto inorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
            if (!node) return;
            self(node->l, out, self);
            out.push_back(node->v);
            self(node->r, out, self);
        };

        std::vector<int> values;
        values.reserve(n);
        for (int i = 1; i <= n; i++) values.push_back(i);
        std::mt19937 rng(seed);
        std::shuffle(values.begin(), values.end(), rng);

        Node *root = nullptr;
        for (int v : values) root = insert(root, v);

        std::vector<int> pre, in;
        preorder(root, pre, preorder);
        inorder(root, in, inorder);
        return {pre, in};
    }
};

// Definition: Encode a tree as P:...;I:... for parity scripts
// Parameters: T: tree to encode
// Returns: string encoding of preorder and inorder
// Errors: returns empty encodings if the tree has empty traversals
static std::string treeToPreInString(const VectorRangeTreeMap &T) {
    auto join = [](const std::vector<int> &vals) {
        std::ostringstream oss;
        for (size_t i = 0; i < vals.size(); i++) {
            if (i) oss << ",";
            oss << vals[i];
        }
        return oss.str();
    };
    std::ostringstream oss;
    oss << "P:" << join(T.original_preorder) << ";I:" << join(T.original_inorder);
    return oss.str();
}

// Definition: Find the smallest k such that FlipDistTree returns true
// Parameters: T_init/T_final: trees; max_k: search cap
// Returns: minimum k, or -1 if not found up to max_k
// Errors: returns -1 on invalid inputs or failures during search
static int FlipDistMinK(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int max_k) {
    if (max_k < 0) return -1;
    const Key128 pair_key = makeKeyPair(T_init, T_final);
    auto &bounds = g_kbounds[pair_key];
    if (bounds.min_success >= 0 && bounds.max_fail >= bounds.min_success) {
        bounds.max_fail = bounds.min_success - 1;
    }

    if (bounds.min_success >= 0 && bounds.max_fail >= bounds.min_success - 1 &&
        bounds.min_success <= max_k) {
        return bounds.min_success;
    }

    int lo = std::max(0, bounds.max_fail + 1);
    // Exact lower bound: every conflicting edge must be removed at least once
    lo = std::max(lo, countConflictEdges(T_init, T_final));
    if (lo > max_k) {
        return -1;
    }

    int hi = -1;
    if (bounds.min_success >= 0 && bounds.min_success <= max_k) {
        hi = bounds.min_success;
    }

    // Find a feasible upper bound progressively instead of probing max_k directly
    if (hi < 0) {
        int probe = lo;
        while (true) {
            if (FlipDistTree(T_init, T_final, probe)) {
                hi = probe;
                if (bounds.min_success < 0 || probe < bounds.min_success) {
                    bounds.min_success = probe;
                }
                break;
            }
            if (probe > bounds.max_fail) bounds.max_fail = probe;
            if (probe >= max_k) {
                return -1;
            }

            long long next_probe = std::max<long long>(probe + 1LL, 2LL * probe + 1LL);
            if (next_probe > max_k) next_probe = max_k;
            probe = static_cast<int>(next_probe);
        }
    }

    if (lo > hi) {
        // Cached bounds can be tighter than the fresh lower bound
        if (bounds.min_success >= 0 && bounds.min_success <= max_k) {
            return bounds.min_success;
        }
        return -1;
    }

    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (FlipDistTree(T_init, T_final, mid)) {
            hi = mid;
            if (bounds.min_success < 0 || mid < bounds.min_success) {
                bounds.min_success = mid;
            }
        } else {
            if (mid > bounds.max_fail) bounds.max_fail = mid;
            lo = mid + 1;
        }
    }

    if (bounds.min_success < 0 || lo < bounds.min_success) {
        bounds.min_success = lo;
    }
    if (lo - 1 > bounds.max_fail) {
        bounds.max_fail = lo - 1;
    }
    return lo;
}

struct CliOptions {
    std::string case_type = "random";
    int n = 0;
    int seed = 0;
    int count = 1;
    int max_k = -1;
    int bfs_cap = 1;
    bool help = false;
};

// Definition: Emit CLI usage help
// Parameters: argv0: program name
// Returns: nothing
// Errors: none
static void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0 << " --case random|comb --n <int> [--seed <int>] [--count <int>]\n"
              << "       [--max-k <int>] [--bfs-cap <int>]\n";
}

// Definition: Parse CLI flags into CliOptions
// Parameters: argc/argv: CLI args; opts: output options
// Returns: true if parsing succeeds, else false
// Errors: prints to stderr on unknown args
static bool parseArgs(int argc, char **argv, CliOptions &opts) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            opts.help = true;
            return true;
        }
        if (arg == "--case" && i + 1 < argc) {
            opts.case_type = argv[++i];
            continue;
        }
        if (arg == "--n" && i + 1 < argc) {
            opts.n = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--seed" && i + 1 < argc) {
            opts.seed = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--count" && i + 1 < argc) {
            opts.count = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--max-k" && i + 1 < argc) {
            opts.max_k = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--bfs-cap" && i + 1 < argc) {
            opts.bfs_cap = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--emit-path" || arg == "--path-ascii") {
            continue; // accepted for compatibility
        }
        std::cerr << "Unknown argument: " << arg << "\n";
        return false;
    }
    return true;
}

// Definition: Execute the FlipDist CLI and emit JSON result lines
// Parameters: argc/argv: CLI args
// Returns: exit code, 0 on success and nonzero on invalid args
// Errors: prints to stderr on invalid args or unsupported cases
static int runCli(int argc, char **argv) {
    CliOptions opts;
    if (!parseArgs(argc, argv, opts)) {
        printUsage(argv[0]);
        return 2;
    }
    if (opts.help || argc == 1) {
        printUsage(argv[0]);
        return opts.help ? 0 : 2;
    }
    if (opts.n <= 0) {
        std::cerr << "Missing or invalid --n\n";
        return 2;
    }
    if (opts.count <= 0) {
        std::cerr << "Invalid --count\n";
        return 2;
    }
    if (opts.case_type != "random" && opts.case_type != "comb") {
        std::cerr << "Invalid --case (use random or comb)\n";
        return 2;
    }
    if (opts.max_k <= 0) {
        opts.max_k = std::max(1, 3 * opts.n + 10);
    }

    for (int i = 0; i < opts.count; i++) {
        int case_seed = opts.case_type == "random" ? (opts.seed + i) : -1;

        std::pair<std::vector<int>, std::vector<int>> pre1_in1;
        std::pair<std::vector<int>, std::vector<int>> pre2_in2;
        if (opts.case_type == "random") {
            pre1_in1 = TreeGenerator::generateRandom(opts.n, case_seed * 2);
            pre2_in2 = TreeGenerator::generateRandom(opts.n, case_seed * 2 + 1);
        } else {
            pre1_in1 = TreeGenerator::generateLeftChain(opts.n);
            pre2_in2 = TreeGenerator::generateRightChain(opts.n);
        }

        VectorRangeTreeMap A, B;
        A.build(pre1_in1.first, pre1_in1.second);
        B.build(pre2_in2.first, pre2_in2.second);

        auto run_one = [&](const VectorRangeTreeMap &T1,
                           const VectorRangeTreeMap &T2,
                           const char *direction) {
            initProfile();
            if (g_profile.enabled) {
                g_profile.start_time = std::chrono::steady_clock::now();
                g_profile.calls_flipdist = 0;
                g_profile.calls_tdi = 0;
                g_profile.calls_tds = 0;
                g_profile.calls_partition = 0;
                g_profile.free_edge_hits = 0;
                g_profile.free_edge_misses = 0;
                g_profile.s_branch_calls = 0;
                g_profile.s_empty_branch_calls = 0;
                g_profile.independent_subsets_initial = 0;
                g_profile.independent_subsets_s = 0;
                g_profile.max_s_size = 0;
                g_profile.max_i_size = 0;
                g_profile.time_flipdist_ms = 0.0;
                g_profile.time_tdi_ms = 0.0;
                g_profile.time_tds_ms = 0.0;
                g_profile.time_partition_ms = 0.0;
            }
            PerformanceTimer timer;
            timer.start();
            resetMemo();
            int dist = FlipDistMinK(T1, T2, opts.max_k);
            double ms = timer.getMilliseconds();

            std::string status = (dist >= 0) ? "ok" : "not_found";
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "{\"case_type\":\"" << opts.case_type << "\""
                      << ",\"n\":" << opts.n
                      << ",\"seed\":" << case_seed
                      << ",\"direction\":\"" << direction << "\""
                      << ",\"distance\":" << dist
                      << ",\"time_ms\":" << ms
                      << ",\"status\":\"" << status << "\""
                      << ",\"tree_a\":\"" << treeToPreInString(T1) << "\""
                      << ",\"tree_b\":\"" << treeToPreInString(T2) << "\""
                      << ",\"max_k\":" << opts.max_k
                      << "}" << std::endl;

            if (g_profile.enabled) {
                std::cout << "[PROFILE] direction=" << direction << " dist=" << dist << "\n"
                          << "  - FlipDistTree calls=" << g_profile.calls_flipdist
                          << " time_ms=" << g_profile.time_flipdist_ms << "\n"
                          << "  - TreeDistI calls=" << g_profile.calls_tdi
                          << " time_ms=" << g_profile.time_tdi_ms << "\n"
                          << "  - TreeDistS calls=" << g_profile.calls_tds
                          << " time_ms=" << g_profile.time_tds_ms << "\n"
                          << "  - partition calls=" << g_profile.calls_partition
                          << " time_ms=" << g_profile.time_partition_ms << "\n"
                          << "  - free_edge hits=" << g_profile.free_edge_hits
                          << " misses=" << g_profile.free_edge_misses << "\n"
                          << "  - S-branch calls=" << g_profile.s_branch_calls
                          << " S-empty calls=" << g_profile.s_empty_branch_calls << "\n"
                          << "  - independent_subsets_initial=" << g_profile.independent_subsets_initial
                          << " max_I_size=" << g_profile.max_i_size << "\n"
                          << "  - max_S_size=" << g_profile.max_s_size << "\n";
            }
        };

        run_one(A, B, "a->b");
        run_one(B, A, "b->a");
    }

    return 0;
}

// Definition: CLI entry point
// Parameters: argc/argv: CLI args
// Returns: process exit code
// Errors: none beyond runCli
int main(int argc, char **argv) {
    return runCli(argc, argv);
}
