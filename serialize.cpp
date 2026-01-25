// Implements hashing, serialisation, and equality helpers for
// `VectorRangeTreeMap`. These routines give the BFS engine stable fingerprints
// (canonical strings, FNV hashes) and lightweight comparisons used throughout
// the solver and comparison harness.

#include "rotation_tree.h"
#include <sstream>
#include <algorithm>
#include <functional>
#include <cstdint>


// 64-bit FNV-1a hash used when serialising trees to stable fingerprints.
static uint64_t fnv1a64(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}

// Convenience wrapper that hashes a tree via its string serialisation.
[[maybe_unused]] static inline uint64_t key64(const VectorRangeTreeMap& T) {
    return fnv1a64(treeToString(T));
}

// Hash combination for integer pairs; shared across the whole project.
size_t PairHash::operator()(const std::pair<int,int>& p) const {
    return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned long long)p.second);
}
bool PairEq::operator()(const std::pair<int,int>& a, const std::pair<int,int>& b) const {
    return a.first==b.first && a.second==b.second;
}

// Structural equality check: compares node sets, ranges, and adjacency.
bool TreesEqual(const VectorRangeTreeMap& A, const VectorRangeTreeMap& B) {
    if (A.original_nodes != B.original_nodes) return false;
    for (int v : A.original_nodes) {
        if (A.getRange(v)     != B.getRange(v))     return false;
        if (A.getLeftChild(v) != B.getLeftChild(v)) return false;
        if (A.getRightChild(v)!= B.getRightChild(v))return false;
        if (A.getParent(v)    != B.getParent(v))    return false;
    }
    return true;
}

// Stable text serialisation that drives canonicalKey/Tree hashers.
std::string treeToString(const VectorRangeTreeMap& T) {
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());
    std::ostringstream oss;
    for (int v : nodes) {
        auto [rs,re] = T.getRange(v);
        int L = T.getLeftChild(v);
        int R = T.getRightChild(v);
        int P = T.getParent(v);
        oss << v << ',' << rs << ',' << re << ',' << L << ',' << R << ',' << P << ';';
    }
    return oss.str();
}

// Canonical traversal-based fingerprint: encodes preorder and inorder sequences.
std::string canonicalTraversalString(const VectorRangeTreeMap& T) {
    auto emit = [](const std::vector<int>& seq) {
        std::ostringstream oss;
        bool first = true;
        for (int v : seq) {
            if (!first) oss << ',';
            first = false;
            oss << v;
        }
        return oss.str();
    };

    std::ostringstream out;
    out << "P:" << emit(T.original_preorder)
        << ";I:" << emit(T.original_inorder);
    return out.str();
}
