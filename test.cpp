#include "A_tree.h"
#include <numeric>
#include <random>
#include <iostream>
#include <chrono>
#include <cassert>
#include <climits>
#include <vector>
#include <queue>
#include <unordered_map>
#include <optional>
#include <string>

struct BFSRun {
    int dist; size_t expanded, enqueued, visited;
    bool timeout, cap_hit; double seconds;
};
BFSRun BFSSearchCapped(const VectorRangeTreeMap& T_s,
                       const VectorRangeTreeMap& T_e,
                       double time_limit_sec = 10.0,
                       size_t visited_cap    = 5'000'000,
                       size_t queue_cap      = 5'000'000);

static int BiBFSSearch(const VectorRangeTreeMap& S, const VectorRangeTreeMap& T,
                       size_t state_cap = 2'000'000)
{
    if (TreesEqual(S, T)) return 0;

    using TV = VectorRangeTreeMap;
    std::queue<TV> qS, qT;
    std::unordered_map<std::string,int> dS, dT;

    auto kS = treeToString(S), kT = treeToString(T);
    dS.emplace(kS, 0); qS.push(S);
    dT.emplace(kT, 0); qT.push(T);

    auto expand = [&](bool fromS)->std::optional<int> {
        auto& q = fromS ? qS : qT;
        auto& mine = fromS ? dS : dT;
        auto& othr = fromS ? dT : dS;

        size_t layer = q.size();
        while (layer--) {
            auto cur = q.front(); q.pop();
            int dcur = mine[treeToString(cur)];
            for (int v : cur.original_nodes) {
                if (cur.getRightChild(v) != TV::NO_CHILD) {
                    auto tmp = cur; tmp.rotateLeft(v);
                    auto k = treeToString(tmp);
                    if (auto it = othr.find(k); it != othr.end()) return dcur + 1 + it->second;
                    if (mine.emplace(k, dcur+1).second) q.push(std::move(tmp));
                }
                if (cur.getLeftChild(v) != TV::NO_CHILD) {
                    auto tmp = cur; tmp.rotateRight(v);
                    auto k = treeToString(tmp);
                    if (auto it = othr.find(k); it != othr.end()) return dcur + 1 + it->second;
                    if (mine.emplace(k, dcur+1).second) q.push(std::move(tmp));
                }
            }
        }
        if (dS.size() + dT.size() > state_cap) return std::optional<int>{INT_MAX};
        return std::nullopt;
    };

    while (!qS.empty() && !qT.empty()) {
        bool fromS = (qS.size() <= qT.size());
        auto hit = expand(fromS);
        if (hit.has_value()) return *hit;
    }
    return INT_MAX;
}


static void buildComb(int n, bool rightComb, std::vector<int>& pre, std::vector<int>& in) {
    in.resize(n); std::iota(in.begin(), in.end(), 1);
    pre.resize(n);
    if (rightComb) std::iota(pre.begin(), pre.end(), 1);
    else for (int i=0;i<n;i++) pre[i] = n - i;
}

struct Shape{ Shape* L{nullptr}; Shape* R{nullptr}; int label{-1}; };
static Shape* buildRandomShape(int n, std::mt19937& rng) {
    if (n<=0) return nullptr; if (n==1) return new Shape();
    std::uniform_int_distribution<int> pick(0, n-1);
    int left_sz = pick(rng), right_sz = n-1-left_sz;
    auto* t = new Shape(); t->L = buildRandomShape(left_sz, rng); t->R = buildRandomShape(right_sz, rng); return t;
}
static void assignInorder(Shape* t, int& next) { if (!t) return; assignInorder(t->L, next); t->label = next++; assignInorder(t->R, next); }
static void dumpInorder(Shape* t, std::vector<int>& in){ if (!t) return; dumpInorder(t->L, in); in.push_back(t->label); dumpInorder(t->R, in); }
static void dumpPreorder(Shape* t, std::vector<int>& pre){ if (!t) return; pre.push_back(t->label); dumpPreorder(t->L, pre); dumpPreorder(t->R, pre); }
static void freeShape(Shape* t){ if (!t) return; freeShape(t->L); freeShape(t->R); delete t; }

static void buildRandomTraversals(int n, std::vector<int>& pre, std::vector<int>& in, std::mt19937& rng) {
    Shape* t = buildRandomShape(n, rng);
    int next = 1; assignInorder(t, next); dumpInorder(t, in); dumpPreorder(t, pre); freeShape(t);
}

// Local LB used for logs
static int LowerBound_MissingTargetEdges(const VectorRangeTreeMap& cur, const VectorRangeTreeMap& tgt) {
    std::unordered_set<std::pair<int,int>,PairHash,PairEq> Et, Ec;
    cur.collectEdges(cur.root, Ec);
    tgt.collectEdges(tgt.root, Et);
    int common = 0; for (auto& e : Ec) if (Et.count(e)) ++common;
    return int(Et.size()) - common;
}

// ---- tests copied from your main() flow
static void test_minimality_small() {
    std::mt19937 rng(12345);
    for (int n=3; n<=12; ++n) {
        for (int t=0; t<10; ++t) {
            std::vector<int> preA,inA,preB,inB;
            buildRandomTraversals(n, preA, inA, rng);
            buildRandomTraversals(n, preB, inB, rng);
            VectorRangeTreeMap A,B; A.build(preA,inA); B.build(preB,inB);
#ifndef NDEBUG
            A.verify(); B.verify();
#endif
            auto t0 = std::chrono::steady_clock::now();
            int d_oracle = BiBFSSearch(A, B, /*cap=*/1'000'000);
            auto t1 = std::chrono::steady_clock::now();
            int d_again  = BiBFSSearch(A, B, /*cap=*/1'000'000);
            auto t2 = std::chrono::steady_clock::now();

            double ms1 = std::chrono::duration<double,std::milli>(t1-t0).count();
            double ms2 = std::chrono::duration<double,std::milli>(t2-t1).count();
            if (d_oracle == INT_MAX) {
                std::cout << "[WARN] cap hit at n=" << n << "\n";
            } else {
                assert(d_oracle == d_again);
                std::cout << "OK n="<<n<<" rand dist="<<d_oracle
                          << " (oracle "<<ms1<<"ms, repeat "<<ms2<<"ms)\n";
            }
        }
        std::vector<int> preL,inL,preR,inR;
        buildComb(n, false, preL, inL);
        buildComb(n, true,  preR, inR);
        VectorRangeTreeMap L,R; L.build(preL,inL); R.build(preR,inR);
#ifndef NDEBUG
        L.verify(); R.verify();
#endif
        int d = BiBFSSearch(L, R, /*cap=*/2'000'000);
        int expected = n-1;
        std::cout << "n="<<n<<" comb->comb dist="<<d<<" expected="<<expected<<"\n";
        if (d != INT_MAX) assert(d == expected);
    }
}

static void capacity_sweep_bibfs_only() {
    const double TIME_LIMIT = 2.0;
    for (int n=8; n<=40; ++n) {
        std::vector<int> preL,inL,preR,inR;
        buildComb(n, false, preL, inL);
        buildComb(n, true,  preR, inR);
        VectorRangeTreeMap L,R; L.build(preL,inL); R.build(preR,inR);
        auto t0 = std::chrono::steady_clock::now();
        int d = BiBFSSearch(L, R, /*cap=*/5'000'000);
        int LB = LowerBound_MissingTargetEdges(L, R);
        auto t1 = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(t1-t0).count();
        std::cout << "[CAP] n="<<n<<" comb dist="<<d
                  << " LB="<<LB<<" UB="<<d
                  << (d!=INT_MAX && d==LB ? " [Optimal: LB=UB]" : "")
                  << " time="<<sec<<"s\n";
        if (sec > TIME_LIMIT || d==INT_MAX) {
            std::cout << "-> BiBFS capacity (worst case) ~ n≈"<<(n-1)<<"\n";
            break;
        }
    }
}

static void capacity_with_my_bfs() {
    for (int n=13; n<=60; ++n) {
        std::vector<int> preL,inL,preR,inR;
        buildComb(n, false, preL, inL);
        buildComb(n, true,  preR, inR);
        VectorRangeTreeMap L,R; L.build(preL,inL); R.build(preR,inR);
        int LB = LowerBound_MissingTargetEdges(L, R);

        auto res = BFSSearchCapped(L, R, /*time_limit_sec=*/3.0,
                                   /*visited_cap=*/5'000'000, /*queue_cap=*/5'000'000);
        std::cout << "[MY-BFS] n="<<n<<" comb dist="<<res.dist
                  << " LB="<<LB<<" UB="<<res.dist
                  << (res.dist!=INT_MAX && res.dist==LB ? " [Optimal: LB=UB]" : "")
                  << " time="<<res.seconds<<"s"
                  << " expanded="<<res.expanded
                  << " visited="<<res.visited
                  << (res.timeout ? " TIMEOUT" : "")
                  << (res.cap_hit ? " CAP" : "")
                  << "\n";
        if (res.timeout || res.cap_hit || res.dist == INT_MAX) break;
    }
}

static void random_with_my_bfs() {
    std::mt19937 rng(42);
    for (int n=13; n<=40; ++n) {
        for (int t=0; t<3; ++t) {
            std::vector<int> preA,inA,preB,inB;
            buildRandomTraversals(n, preA, inA, rng);
            buildRandomTraversals(n, preB, inB, rng);
            VectorRangeTreeMap A,B; A.build(preA,inA); B.build(preB,inB);
            int LB = LowerBound_MissingTargetEdges(A, B);

            auto res = BFSSearchCapped(A, B, 2.0, 5'000'000, 5'000'000);
            std::cout << "   n="<<n<<" rand dist="<<res.dist
                      << " LB="<<LB<<" UB="<<res.dist
                      << (res.dist!=INT_MAX && res.dist==LB ? " [Optimal: LB=UB]" : "")
                      << " time="<<res.seconds<<"s "
                      << (res.timeout ? "TIMEOUT " : "")
                      << (res.cap_hit ? "CAP " : "")
                      << "\n";
            if (res.timeout || res.cap_hit) return;
        }
    }
}

void runAllTests() {
    test_minimality_small();
    capacity_sweep_bibfs_only();
    capacity_with_my_bfs();
    random_with_my_bfs();
}

int main()
{
    runAllTests();
    return 0;
}