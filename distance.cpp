#include "A_tree.h"
#include <unordered_map>
#include <climits>
#include <cassert>

// memo + recursion + remove

static std::unordered_map<std::string,int> memo;

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

int FindRotationDistance(const VectorRangeTreeMap& T_init,
                         const VectorRangeTreeMap& T_final)
{
    memo.clear();
    return Dist(T_init, T_final);
}
