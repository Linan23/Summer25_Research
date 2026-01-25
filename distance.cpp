// Implements higher-level distance helpers that combine free-edge
// decompositions, memoisation, and fallbacks to the BFS engine. These helpers
// back the heuristics used by the optimised solver when it estimates lower
// bounds or short-circuits trivially solvable subproblems.

#include "rotation_tree.h"
#include <unordered_map>
#include <climits>
#include <cassert>

// Convenience wrapper used by heuristics to count active nodes.
static inline size_t node_count(const VectorRangeTreeMap& T) {
    return T.original_nodes.size();
}

// memo + recursion + remove


// Standard 64-bit FNV-1a hash reused by the memoisation tables.
static uint64_t fnv1a64(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}

// 64-bit FNV 
// Hashes a tree via its canonical string; used by the Dist memo tables.
static uint64_t hashTree64(const VectorRangeTreeMap& T) {
    return fnv1a64(treeToString(T));
}

struct Key2 {
    uint64_t a, b;
    bool operator==(const Key2& o) const noexcept { return a == o.a && b == o.b; }
};

struct Key2H {
    size_t operator()(const Key2& k) const noexcept {
        return size_t(k.a ^ (k.b * 0x9e3779b97f4a7c15ULL));
    }
};

static std::unordered_map<Key2,int,Key2H> M;


static std::unordered_map<std::string,int> memo;

/*
int Dist(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te) {
    if (TreesEqual(Ts, Te)) return 0;

    std::string key = treeToString(Ts) + "|" + treeToString(Te);
    if (auto it = memo.find(key); it != memo.end()) return it->second;

    int u; bool left;
    if (hasFreeEdge(Ts, Te, u, left)) {
        VectorRangeTreeMap Tbar = Ts;
        if (left) Tbar.rotateLeft(u); else Tbar.rotateRight(u);

        int v = left ? Tbar.getRightChild(u) : Tbar.getLeftChild(u);
        auto pr = Tbar.getRange(u);
        auto cr = Tbar.getRange(v);

        auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(Tbar, pr, cr);
        auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(Te,   pr, cr);

        int result;
        if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes) {
            result = BFSSearch(Ts, Te);
        } else {
            result = 1 + Dist(A1, B1) + Dist(A2, B2);
            if (result == 1) {
                auto Tcheck = Ts;
                if (left) Tcheck.rotateLeft(u); else Tcheck.rotateRight(u);
                if (!TreesEqual(Tcheck, Te)) result = BFSSearch(Ts, Te);
            }
        }
        memo[key] = result; return result;
    } else {
        int result = BFSSearch(Ts, Te);
        memo[key] = result; return result;
    }
}
*/

// Depth-limited distance calculation that falls back to full BFS once the
[[maybe_unused]] static int DistRec(const VectorRangeTreeMap& T_s,
                                    const VectorRangeTreeMap& T_e,
                                    int budget) {
    if (TreesEqual(T_s, T_e)) return 0;

    Key2 key{hashTree64(T_s), hashTree64(T_e)};
    if (auto it = M.find(key); it != M.end()) return it->second;

    // Budget exhausted? exact BFS
    if (budget <= 0) {
        int d = BFSSearch(T_s, T_e);
        M[key] = d; return d;
    }

    int u; bool left;
    if (hasFreeEdge(T_s, T_e, u, left)) {
        // removeFreeEdge will call Dist again; pass reduced budget
        int d = removeFreeEdge(T_s, T_e);
        M[key] = d; return d;
    }

    int d_bfs = BFSSearch(T_s, T_e);
    M[key] = d_bfs; return d_bfs;
}


// Recursive wrapper implementing memoised distance computation with a fallback
// to BFS when the recursion gets too deep.
static int DistImpl(const VectorRangeTreeMap& Ts,
                    const VectorRangeTreeMap& Te,
                    int depth, int depth_cap)
{
    if (TreesEqual(Ts, Te)) return 0;

    std::string key = treeToString(Ts) + "|" + treeToString(Te);
    if (auto it = memo.find(key); it != memo.end()) return it->second;

    // depth cap: fall back to brute BFS if recursion is getting silly
    if (depth > depth_cap) {
        int ans = BFSSearch(Ts, Te);
        memo.emplace(std::move(key), ans);
        return ans;
    }

    // mark in-flight to break any accidental cycles
    memo.emplace(key, INT_MAX);

    int best = INT_MAX;
    int u; bool left;
    if (hasFreeEdge(Ts, Te, u, left)) {
        int try_free = removeFreeEdge(Ts, Te);         // no BFSSearch inside
        if (try_free != INT_MAX) best = try_free;
    }
    if (best == INT_MAX) {
        best = BFSSearch(Ts, Te);                      // brute once
    }

    memo[key] = best;
    return best;
}

// Public entry point for the distance computation with heuristic pruning.
int Dist(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te) {
    const int depth_cap = (int)(2 * std::max(node_count(Ts), node_count(Te)) + 4);
    return DistImpl(Ts, Te, /*depth=*/0, depth_cap);
}



/*
// Applies a “free edge” rotation when it provably reduces the problem into two
// smaller subproblems; otherwise returns INT_MAX to signal fallback to BFS.
int removeFreeEdge(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te) {
    int u; bool left;
    assert(hasFreeEdge(Ts, Te, u, left));

    VectorRangeTreeMap Tbar = Ts;
    if (left) Tbar.rotateLeft(u); else Tbar.rotateRight(u);

    int v = left ? Tbar.getRightChild(u) : Tbar.getLeftChild(u);
    auto pr = Tbar.getRange(u);
    auto cr = Tbar.getRange(v);

    auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(Tbar, pr, cr);
    auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(Te,   pr, cr);

    if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes)
        return BFSSearch(Ts, Te);

    return 1 + Dist(A1, B1) + Dist(A2, B2);
}
*/

int removeFreeEdge(const VectorRangeTreeMap& Ts, const VectorRangeTreeMap& Te) {
    int u; bool left;
    assert(hasFreeEdge(Ts, Te, u, left));

    VectorRangeTreeMap Tbar = Ts;
    if (left) Tbar.rotateLeft(u); else Tbar.rotateRight(u);

    // sanity: the rotation must have created some target edge
    {
        std::unordered_set<std::pair<int,int>,PairHash,PairEq> Et, Eb;
        Te .collectEdges(Te.root, Et);
        Tbar.collectEdges(Tbar.root, Eb);
        bool created_common = false;
        for (auto const& e : Eb) {
            if (!Ts.isOriginal(e.first) || !Ts.isOriginal(e.second)) continue;
            // if this edge wasn't in Ts but is in Te, we created a target edge
            std::unordered_set<std::pair<int,int>,PairHash,PairEq> Es;
            Ts.collectEdges(Ts.root, Es);
            if (!Es.count(e) && Et.count(e)) { created_common = true; break; }
        }
        if (!created_common) return INT_MAX; // not a valid free-edge shortcut
    }

    int v  = left ? Tbar.getRightChild(u) : Tbar.getLeftChild(u);
    auto pr = Tbar.getRange(u);
    auto cr = Tbar.getRange(v);

    auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(Tbar, pr, cr);
    auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(Te,   pr, cr);

    // partitions must match AND strictly reduce problem size
    const size_t N  = node_count(Ts);
    const size_t n1 = node_count(A1), n2 = node_count(A2);
    if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes)
        return INT_MAX;
    if (n1 == 0 || n2 == 0) return INT_MAX;
    if (n1 >= N || n2 >= N) return INT_MAX;

    // 1 rotation for the created edge + subproblems
    int d1 = Dist(A1, B1); if (d1 == INT_MAX) return INT_MAX;
    int d2 = Dist(A2, B2); if (d2 == INT_MAX) return INT_MAX;
    return 1 + d1 + d2;
}

// High-level API used by the CLI/tests: clears memo state and runs Dist.
int FindRotationDistance(const VectorRangeTreeMap& T_init,
                         const VectorRangeTreeMap& T_final)
{
    memo.clear();
    return Dist(T_init, T_final);
}
