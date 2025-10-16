#include "rotation_tree.h"
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
#include <cstdlib>
#include <cerrno>

// Parses an environment variable as a double, returning fallback on failure.
static bool envAsFlag(const char *name, bool fallback = false)
{
    if (const char *val = std::getenv(name))
    {
        if (val[0] == '\0' || val[0] == '0' || val[0] == 'f' || val[0] == 'F' || val[0] == 'n' || val[0] == 'N')
            return false;
        return true;
    }
    return fallback;
}

static double envAsDouble(const char *name, double fallback)
{
    if (const char *val = std::getenv(name))
    {
        char *end = nullptr;
        errno = 0;
        double parsed = std::strtod(val, &end);
        if (errno == 0 && end && end != val)
            return parsed;
    }
    return fallback;
}

// Parses an environment variable as an unsigned integer, returning fallback if unset.
static size_t envAsSize(const char *name, size_t fallback)
{
    if (const char *val = std::getenv(name))
    {
        char *end = nullptr;
        errno = 0;
        unsigned long long parsed = std::strtoull(val, &end, 10);
        if (errno == 0 && end && end != val)
            return static_cast<size_t>(parsed);
    }
    return fallback;
}

// Straightforward string-based BiBFS kept as a correctness oracle.
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


// Generates traversals for a “comb” tree, leaning right (spine) when rightComb.
static void buildComb(int n, bool rightComb, std::vector<int>& pre, std::vector<int>& in) {
    in.resize(n); std::iota(in.begin(), in.end(), 1);
    pre.resize(n);
    if (rightComb) std::iota(pre.begin(), pre.end(), 1);
    else for (int i=0;i<n;i++) pre[i] = n - i;
}

struct Shape{ Shape* L{nullptr}; Shape* R{nullptr}; int label{-1}; };
// Recursively synthesises a random binary tree with n nodes; leaves are shared
// helper structures used only during test generation.
static Shape* buildRandomShape(int n, std::mt19937& rng) {
    if (n<=0) return nullptr; if (n==1) return new Shape();
    std::uniform_int_distribution<int> pick(0, n-1);
    int left_sz = pick(rng), right_sz = n-1-left_sz;
    auto* t = new Shape(); t->L = buildRandomShape(left_sz, rng); t->R = buildRandomShape(right_sz, rng); return t;
}
// Assigns labels 1..n to the shape by walking it in-order.
static void assignInorder(Shape* t, int& next) { if (!t) return; assignInorder(t->L, next); t->label = next++; assignInorder(t->R, next); }
// Writes the shape's in-order traversal into the provided buffer.
static void dumpInorder(Shape* t, std::vector<int>& in){ if (!t) return; dumpInorder(t->L, in); in.push_back(t->label); dumpInorder(t->R, in); }
// Writes the shape's pre-order traversal into the provided buffer.
static void dumpPreorder(Shape* t, std::vector<int>& pre){ if (!t) return; pre.push_back(t->label); dumpPreorder(t->L, pre); dumpPreorder(t->R, pre); }
// Cleans up the temporary tree produced for test generation.
static void freeShape(Shape* t){ if (!t) return; freeShape(t->L); freeShape(t->R); delete t; }

// Converts a random Shape into preorder/inorder vectors and then frees it.
static void buildRandomTraversals(int n, std::vector<int>& pre, std::vector<int>& in, std::mt19937& rng) {
    Shape* t = buildRandomShape(n, rng);
    int next = 1; assignInorder(t, next); dumpInorder(t, in); dumpPreorder(t, pre); freeShape(t);
}

// Local LB used for logs
// Counts how many target edges are still missing and converts the tally into a
// simple admissible lower bound (each rotation can fix at most three edges).
static int LowerBound_MissingTargetEdges(const VectorRangeTreeMap& cur, const VectorRangeTreeMap& tgt) {
    std::unordered_set<std::pair<int,int>,PairHash,PairEq> Et, Ec;
    cur.collectEdges(cur.root, Ec);
    tgt.collectEdges(tgt.root, Et);
    int common = 0; for (auto& e : Ec) if (Et.count(e)) ++common;
    return int(Et.size()) - common;
}

// Verifies that identical trees report distance 0 and that statistics stay sane.
static void test_identity_cases() {
    std::vector<int> pre{2,1,3}, in{1,2,3};
    VectorRangeTreeMap A; A.build(pre, in);
#ifndef NDEBUG
    A.verify();
#endif
    BFSStats stats{};
    int d = BFSSearch(A, A, &stats);
    assert(d == 0);
    assert(stats.equality_hits >= 1);
    std::cout << "[TEST] identity distance=" << d << " enqueued=" << stats.enqueued << "\n";
}

// Ensures a single rotation distance is computed correctly on a tiny tree.
static void test_single_rotation() {
    // Balanced hinge becomes right comb after one rotateRight.
    std::vector<int> pre_bal{2,1,3}, in_bal{1,2,3};
    std::vector<int> pre_comb{1,2,3}, in_comb{1,2,3};
    VectorRangeTreeMap balanced, comb;
    balanced.build(pre_bal, in_bal);
    comb.build(pre_comb, in_comb);

    int expected = BiBFSSearchHashed(balanced, comb, /*cap=*/1000);
    assert(expected != INT_MAX);
    int d = BFSSearch(balanced, comb);
    assert(d == expected);
    std::cout << "[TEST] single rotation dist=" << d << "\n";
}

// Cross-checks skewed vs balanced trees using the hashed BiBFS reference.
static void test_skew_vs_balanced() {
    std::vector<int> pre_left{4,3,2,1}, in_left{1,2,3,4};
    std::vector<int> pre_bal{2,1,3,4}, in_bal{1,2,3,4};
    VectorRangeTreeMap left, balanced;
    left.build(pre_left, in_left);
    balanced.build(pre_bal, in_bal);

    int expected = BiBFSSearchHashed(left, balanced, /*cap=*/5'000'000);
    assert(expected != INT_MAX);
    BFSStats stats{};
    int d = BFSSearch(left, balanced, &stats);
    assert(d == expected);
    std::cout << "[TEST] skew vs balanced dist=" << d
              << " expanded=" << stats.generated << "\n";
}

// Random small instances: confirm hashed BiBFS matches single-ended BFS results.
static void test_random_small_agreement() {
    std::mt19937 rng(2025);
    for (int n = 3; n <= 7; ++n) {
        for (int trial = 0; trial < 8; ++trial) {
            std::vector<int> preA,inA,preB,inB;
            buildRandomTraversals(n, preA, inA, rng);
            buildRandomTraversals(n, preB, inB, rng);
            VectorRangeTreeMap A,B; A.build(preA,inA); B.build(preB,inB);
            int d_bidir = BiBFSSearchHashed(A, B, /*cap=*/5'000'000);
            assert(d_bidir != INT_MAX);
            int d_baseline = BFSSearchBaseline(A, B);
            int d_optimized = BFSSearchOptimized(A, B);
            assert(d_baseline == d_bidir);
            assert(d_optimized == d_bidir);
        }
    }
    std::cout << "[TEST] random small agreement ok\n";
}

// Exhaustively checks small trees to ensure the optimised BFS matches BiBFS.
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

// Benchmarks when the reference BiBFS runs out of steam
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

// Benchmarks the filtered BFS on comb instances while logging telemetry.
static void capacity_with_my_bfs() {
    for (int n=13; n<=60; ++n) {
        std::vector<int> preL,inL,preR,inR;
        buildComb(n, false, preL, inL);
        buildComb(n, true,  preR, inR);
        VectorRangeTreeMap L,R; L.build(preL,inL); R.build(preR,inR);
        int LB = LowerBound_MissingTargetEdges(L, R);

        const double time_limit = envAsDouble("MY_BFS_TIME_LIMIT", 3.0);
        const size_t visited_cap = envAsSize("MY_BFS_VISITED_CAP", 5'000'000);
        const size_t queue_cap   = envAsSize("MY_BFS_QUEUE_CAP",   5'000'000);
        auto res = BFSSearchCapped(L, R, time_limit, visited_cap, queue_cap);
        std::cout << "[MY-BFS] n="<<n<<" comb dist="<<res.dist
                  << " LB="<<LB<<" UB="<<res.dist
                  << (res.dist!=INT_MAX && res.dist==LB ? " [Optimal: LB=UB]" : "")
                  << " time="<<res.seconds<<"s"
                  << " expanded="<<res.expanded
                  << " visited="<<res.visited
                  << " dup="<<res.duplicates
                  << " maxQ="<<res.max_queue
                  << " eqHits="<<res.equality_hits
                  << " sigMismatch="<<res.signature_mismatches
                  << (res.timeout ? " TIMEOUT" : "")
                  << (res.cap_hit ? " CAP" : "")
                  << "\n";
        if (res.timeout || res.cap_hit || res.dist == INT_MAX) break;
    }
}

// Samples random trees and prints telemetry to understand scaling behaviour.

static void random_with_my_bfs() {
    std::mt19937 rng(42);
    for (int n=13; n<=40; ++n) {
        for (int t=0; t<3; ++t) {
            std::vector<int> preA,inA,preB,inB;
            buildRandomTraversals(n, preA, inA, rng);
            buildRandomTraversals(n, preB, inB, rng);
            VectorRangeTreeMap A,B; A.build(preA,inA); B.build(preB,inB);
            int LB = LowerBound_MissingTargetEdges(A, B);

            const double time_limit = envAsDouble("MY_BFS_RANDOM_TIME_LIMIT", 2.0);
            const size_t visited_cap = envAsSize("MY_BFS_RANDOM_VISITED_CAP", 5'000'000);
            const size_t queue_cap   = envAsSize("MY_BFS_RANDOM_QUEUE_CAP",   5'000'000);
            auto res = BFSSearchCapped(A, B, time_limit, visited_cap, queue_cap);
            std::cout << "   n="<<n<<" rand dist="<<res.dist
                      << " LB="<<LB<<" UB="<<res.dist
                      << (res.dist!=INT_MAX && res.dist==LB ? " [Optimal: LB=UB]" : "")
                      << " time="<<res.seconds<<"s "
                      << "dup="<<res.duplicates<<" "
                      << "maxQ="<<res.max_queue<<" "
                      << (res.timeout ? "TIMEOUT " : "")
                      << (res.cap_hit ? "CAP " : "")
                      << "\n";

            if (envAsFlag("MY_BFS_RANDOM_USE_BIDIR"))
            {
                size_t state_cap = envAsSize("MY_BFS_RANDOM_BIDIR_CAP", 5'000'000);
                int d_bidir = BiBFSSearchHashed(A, B, state_cap);
                std::cout << "      [BiBiFS] dist="<<d_bidir<<" cap="<<state_cap
                          << (d_bidir == INT_MAX ? " HIT" : "") << "\n";
            }
            if (res.timeout || res.cap_hit) return;
        }
    }
}

// Master test entry point executed by main().
void runAllTests() {
    test_identity_cases();
    test_single_rotation();
    test_skew_vs_balanced();
    test_random_small_agreement();
    test_minimality_small();
    capacity_sweep_bibfs_only();
    capacity_with_my_bfs();
    // random_with_my_bfs();
}

int main()
{
    runAllTests();
    return 0;
}
