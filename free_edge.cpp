// Contains utilities for detecting and resolving “free edges” (edges that can
// be rotated independently) during distance computations. By splitting the
// tree along these edges we reduce the search space the BFS needs to explore
// and cache subproblem distances.

#include "rotation_tree.h"
#include <unordered_set>
#include <functional>

// RP hashing/equality
bool RP::operator==(RP const& o) const { return ps==o.ps && pe==o.pe && cs==o.cs && ce==o.ce; }
size_t RPH::operator()(RP const& r) const {
    uint64_t h1 = ((uint64_t)r.ps << 32) ^ uint32_t(r.pe);
    uint64_t h2 = ((uint64_t)r.cs << 32) ^ uint32_t(r.ce);
    return std::hash<uint64_t>()(h1 ^ (h2 * 0x9e3779b97f4a7c15ULL));
}

// Collects the range-based edge signature for every edge in the target tree.
std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap& T) {
    std::unordered_set<RP,RPH> s;
    std::function<void(int)> dfs = [&](int v){
        if (v < 0 || !T.isOriginal(v)) return;
        auto pr = T.getRange(v);
        for (int c : {T.getLeftChild(v), T.getRightChild(v)}) {
            if (c >= 0 && T.isOriginal(c)) {
                auto cr = T.getRange(c);
                s.insert(RP{pr.first, pr.second, cr.first, cr.second});
                dfs(c);
            }
        }
    };
    dfs(T.root);
    return s;
}

// Returns true if the given range-based edge exists in tree.
bool hasEdgeByRange(const VectorRangeTreeMap& tree, const RP& e) {
    std::function<bool(int)> dfs = [&](int v)->bool{
        if (v < 0 || !tree.isOriginal(v)) return false;
        auto pr = tree.getRange(v);
        if (pr.first==e.ps && pr.second==e.pe) {
            for (int c : {tree.getLeftChild(v), tree.getRightChild(v)}) {
                if (c >= 0 && tree.isOriginal(c)) {
                    auto cr = tree.getRange(c);
                    if (cr.first==e.cs && cr.second==e.ce) return true;
                }
            }
        }
        return dfs(tree.getLeftChild(v)) || dfs(tree.getRightChild(v));
    };
    return dfs(tree.root);
}

// internal
static void collectEdgeSet(const VectorRangeTreeMap& T,
    std::unordered_set<std::pair<int,int>,PairHash,PairEq>& S)
{
    std::function<void(int)> dfs = [&](int v){
        if (v < 0 || !T.isOriginal(v)) return;
        int L = T.getLeftChild(v), R = T.getRightChild(v);
        if (L >= 0 && T.isOriginal(L)) { S.insert({v,L}); dfs(L); }
        if (R >= 0 && T.isOriginal(R)) { S.insert({v,R}); dfs(R); }
    };
    dfs(T.root);
}

// Searches for a rotation in cur that creates an edge already present in tgt.
bool hasFreeEdge(const VectorRangeTreeMap& cur,
                 const VectorRangeTreeMap& tgt,
                 int& out_v, bool& out_leftRotation)
{
    std::unordered_set<std::pair<int,int>,PairHash,PairEq> tgtEdges, curEdges;
    collectEdgeSet(tgt, tgtEdges);
    collectEdgeSet(cur, curEdges);

    auto creates_target_edge = [&](int v, bool left)->bool{
        auto t = cur;
        if (left) t.rotateLeft(v); else t.rotateRight(v);
        std::unordered_set<std::pair<int,int>,PairHash,PairEq> tEdges;
        collectEdgeSet(t, tEdges);
        for (auto const& e : tEdges) {
            if (!curEdges.count(e) && tgtEdges.count(e)) return true;
        }
        return false;
    };

    for (int v : cur.original_nodes) {
        if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
            if (creates_target_edge(v, true))  { out_v = v; out_leftRotation = true;  return true; }
        }
        if (cur.getLeftChild(v)  != VectorRangeTreeMap::NO_CHILD) {
            if (creates_target_edge(v, false)) { out_v = v; out_leftRotation = false; return true; }
        }
    }
    return false;
}
