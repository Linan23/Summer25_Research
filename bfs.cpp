#include "A_tree.h"
#include <queue>
#include <unordered_set>
#include <climits>
#include <chrono>

int BFSSearch(const VectorRangeTreeMap& T_s, const VectorRangeTreeMap& T_e) {
    if (TreesEqual(T_s, T_e)) return 0;

    std::queue<std::pair<VectorRangeTreeMap,int>> Q;
    std::unordered_set<std::string> visited;

    Q.push({T_s, 0});
    visited.insert(treeToString(T_s));

    while (!Q.empty()) {
        auto [cur, d] = Q.front(); Q.pop();
        int nd = d + 1;

        for (int v : cur.original_nodes) {
            if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur; tmp.rotateLeft(v);
                if (TreesEqual(tmp, T_e)) return nd;

                int u2; bool l2;
                if (hasFreeEdge(tmp, T_e, u2, l2)) return nd + removeFreeEdge(tmp, T_e);

                auto key = treeToString(tmp);
                if (!visited.count(key)) { visited.insert(key); Q.push({std::move(tmp), nd}); }
            }
            if (cur.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur; tmp.rotateRight(v);
                if (TreesEqual(tmp, T_e)) return nd;

                int u2; bool l2;
                if (hasFreeEdge(tmp, T_e, u2, l2)) return nd + removeFreeEdge(tmp, T_e);

                auto key = treeToString(tmp);
                if (!visited.count(key)) { visited.insert(key); Q.push({std::move(tmp), nd}); }
            }
        }
    }
    return INT_MAX;
}


struct BFSRun {
    int dist; size_t expanded, enqueued, visited;
    bool timeout, cap_hit; double seconds;
};

static uint64_t fnv1a64(const std::string& s){
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}
static inline uint64_t key64(const VectorRangeTreeMap& T) { return fnv1a64(treeToString(T)); }

BFSRun BFSSearchCapped(const VectorRangeTreeMap& T_s,
                              const VectorRangeTreeMap& T_e,
                              double time_limit_sec = 10.0,
                              size_t visited_cap = 5'000'000,
                              size_t queue_cap   = 5'000'000)
{
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();

    std::queue<std::pair<VectorRangeTreeMap,int>> Q;
    std::unordered_set<uint64_t> V;
    V.reserve(std::min(visited_cap, size_t(2'000'000)));

    size_t expanded = 0, enqueued = 0;

    auto push_state = [&](VectorRangeTreeMap&& X, int d)->bool{
        if (V.size() >= visited_cap || Q.size() >= queue_cap) return false;
        uint64_t k = key64(X);
        if (V.emplace(k).second) { Q.push({std::move(X), d}); ++enqueued; return true; }
        return true;
    };

    if (!push_state(VectorRangeTreeMap(T_s), 0)) {
        return {INT_MAX, 0, 0, 0, false, true, 0.0};
    }

    while (!Q.empty()) {
        auto now = Clock::now();
        double sec = std::chrono::duration<double>(now - t0).count();
        if (sec > time_limit_sec) return {INT_MAX, expanded, enqueued, V.size(), true, false, sec};

        auto [cur, d] = Q.front(); Q.pop();
        ++expanded; int nextD = d + 1;

        for (int v : cur.original_nodes) {
            if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur; tmp.rotateLeft(v);
                if (TreesEqual(tmp, T_e))
                    return {nextD, expanded, enqueued+1, V.size(), false, false,
                            std::chrono::duration<double>(Clock::now() - t0).count()};
                if (!push_state(std::move(tmp), nextD))
                    return {INT_MAX, expanded, enqueued, V.size(), false, true,
                            std::chrono::duration<double>(Clock::now() - t0).count()};
            }
            if (cur.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur; tmp.rotateRight(v);
                if (TreesEqual(tmp, T_e))
                    return {nextD, expanded, enqueued+1, V.size(), false, false,
                            std::chrono::duration<double>(Clock::now() - t0).count()};
                if (!push_state(std::move(tmp), nextD))
                    return {INT_MAX, expanded, enqueued, V.size(), false, true,
                            std::chrono::duration<double>(Clock::now() - t0).count()};
            }
        }
    }
    return {INT_MAX, expanded, enqueued, V.size(), false, false,
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()};
}
