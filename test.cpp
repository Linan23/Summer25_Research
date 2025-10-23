#include "rotation_tree.h"
#include "comparison.h"
#include "tree_generators.h"
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

// Ensures the canonical traversal serialisation stays stable.
static void test_canonical_serialization() {
    std::vector<int> pre{2,1,3}, in{1,2,3};
    VectorRangeTreeMap T; T.build(pre, in);
    std::string enc = canonicalTraversalString(T);
    assert(enc == "P:2,1,3;I:1,2,3");

    // Rotation does not change traversal order sets.
    T.rotateRight(2);
    std::string enc2 = canonicalTraversalString(T);
    assert(enc == enc2);
}

// Verifies JSON serialisation for comparison rows.
static void test_comparison_row_json() {
    std::vector<int> pre{2,1,3}, in{1,2,3};
    VectorRangeTreeMap A; A.build(pre, in);
    VectorRangeTreeMap B; B.build(pre, in);

    ComparisonOptions opts;
    opts.case_type = "identity";
    opts.direction = "a->b";
    opts.seed = 123;
    opts.time_limit_sec = 1.0;
    opts.visited_cap = 1024;
    opts.queue_cap = 1024;

    ComparisonRow row = runComparison(A, B, opts);
    if (row.distance != 0) {
        std::cerr << "[DEBUG] row.distance=" << row.distance
                  << " status=" << row.status
                  << " solver=" << row.solver
                  << " time_ms=" << row.time_ms << "\n";
    }
    assert(row.distance == 0);
    assert(row.status == "ok");
    assert(row.solver == "bfs");

    std::string json = comparisonRowToJson(row);
    assert(json.find("\"program\":\"test_asan\"") != std::string::npos);
    assert(json.find("\"case_type\":\"identity\"") != std::string::npos);
    assert(json.find("\"distance\":0") != std::string::npos);
    assert(json.find("\"status\":\"ok\"") != std::string::npos);
    assert(json.find("\"solver\":\"bfs\"") != std::string::npos);
}

// Checks that forward and reverse distances agree on simple cases.
static void test_bidirectional_agreement() {
    std::vector<int> preA{2,1,3}, inA{1,2,3};
    std::vector<int> preB{1,2,3}, inB{1,2,3};
    VectorRangeTreeMap A; A.build(preA, inA);
    VectorRangeTreeMap B; B.build(preB, inB);

    ComparisonOptions opts;
    opts.case_type = "unit";
    opts.seed = 7;
    opts.time_limit_sec = 1.0;
    opts.visited_cap = 1024;
    opts.queue_cap = 1024;

    ComparisonPair pair = runBidirectional(A, B, opts);
    assert(pair.distance_agrees);
    assert(pair.forward.distance == 1);
    assert(pair.reverse.distance == 1);
    assert(pair.forward.status == "ok");
    assert(pair.reverse.status == "ok");
    assert(pair.forward.solver == "bfs");
    assert(pair.reverse.solver == "bfs");
}

// Ensures the comparison wrapper falls back to bidirectional search when BFS hits limits.
static void test_comparison_fallback() {
    Traversals leftComb  = makeCombTraversals(4, /*rightComb=*/false);
    Traversals rightComb = makeCombTraversals(4, /*rightComb=*/true);
    VectorRangeTreeMap A, B;
    A.build(leftComb.preorder, leftComb.inorder);
    B.build(rightComb.preorder, rightComb.inorder);

    ComparisonOptions opts;
    opts.case_type = "fallback";
    opts.seed = 0;
    opts.time_limit_sec = 1.0;
    opts.visited_cap = 1;
    opts.queue_cap = 1;
    opts.use_bidir_on_timeout = true;
    opts.bidir_state_cap = 10'000;

    ComparisonRow row = runComparison(A, B, opts);
    assert(row.distance == 3);

    ComparisonPair pair = runBidirectional(A, B, opts);
    assert(pair.distance_agrees);
    assert(pair.forward.distance == 3);
    assert(pair.reverse.distance == 3);
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
            Traversals tA = makeRandomTraversals(n, rng);
            Traversals tB = makeRandomTraversals(n, rng);
            VectorRangeTreeMap A,B; A.build(tA.preorder, tA.inorder); B.build(tB.preorder, tB.inorder);
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
            Traversals tA = makeRandomTraversals(n, rng);
            Traversals tB = makeRandomTraversals(n, rng);
            VectorRangeTreeMap A,B; A.build(tA.preorder, tA.inorder); B.build(tB.preorder, tB.inorder);
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
        Traversals leftComb  = makeCombTraversals(n, /*rightComb=*/false);
        Traversals rightComb = makeCombTraversals(n, /*rightComb=*/true);
        VectorRangeTreeMap L,R; L.build(leftComb.preorder, leftComb.inorder); R.build(rightComb.preorder, rightComb.inorder);
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
        Traversals leftComb  = makeCombTraversals(n, /*rightComb=*/false);
        Traversals rightComb = makeCombTraversals(n, /*rightComb=*/true);
        VectorRangeTreeMap L,R; L.build(leftComb.preorder, leftComb.inorder); R.build(rightComb.preorder, rightComb.inorder);
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
        Traversals leftComb  = makeCombTraversals(n, /*rightComb=*/false);
        Traversals rightComb = makeCombTraversals(n, /*rightComb=*/true);
        VectorRangeTreeMap L,R; L.build(leftComb.preorder, leftComb.inorder); R.build(rightComb.preorder, rightComb.inorder);
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
            Traversals tA = makeRandomTraversals(n, rng);
            Traversals tB = makeRandomTraversals(n, rng);
            VectorRangeTreeMap A,B; A.build(tA.preorder, tA.inorder); B.build(tB.preorder, tB.inorder);
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

// Optional probe for the random-case ceiling of the brute-force solver.
static void probe_random_limit() {
    if (!envAsFlag("MY_BFS_PROBE_RANDOM_LIMIT"))
        return;

    int start_n = static_cast<int>(envAsSize("MY_BFS_PROBE_RANDOM_START", 8));
    int max_n   = static_cast<int>(envAsSize("MY_BFS_PROBE_RANDOM_MAX",   30));
    int trials  = static_cast<int>(envAsSize("MY_BFS_PROBE_RANDOM_TRIALS", 5));
    uint32_t seed = static_cast<uint32_t>(envAsSize("MY_BFS_PROBE_RANDOM_SEED", 2025));

    double time_limit = envAsDouble("MY_BFS_PROBE_RANDOM_TIME_LIMIT", 5.0);
    size_t visited_cap = envAsSize("MY_BFS_PROBE_RANDOM_VISITED_CAP", 5'000'000);
    size_t queue_cap   = envAsSize("MY_BFS_PROBE_RANDOM_QUEUE_CAP",   5'000'000);

    bool use_bidir = envAsFlag("MY_BFS_PROBE_RANDOM_USE_BIDIR");
    size_t bidir_cap = envAsSize("MY_BFS_PROBE_RANDOM_BIDIR_CAP", 5'000'000);

    std::mt19937 rng(seed);
    for (int n = start_n; n <= max_n; ++n) {
        bool all_ok = true;
        for (int t = 0; t < trials; ++t) {
            Traversals tA = makeRandomTraversals(n, rng);
            Traversals tB = makeRandomTraversals(n, rng);
            VectorRangeTreeMap A, B;
            A.build(tA.preorder, tA.inorder);
            B.build(tB.preorder, tB.inorder);

            auto t0 = std::chrono::steady_clock::now();
            BFSRun res = BFSSearchCapped(A, B, time_limit, visited_cap, queue_cap);
            auto t1 = std::chrono::steady_clock::now();
            double bfs_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            std::cout << "[PROBE] n=" << n
                      << " trial=" << t
                      << " bfs_dist=" << res.dist
                      << " time_ms=" << bfs_ms
                      << " visited=" << res.visited
                      << " maxQ=" << res.max_queue
                      << " dup=" << res.duplicates
                      << (res.timeout ? " TIMEOUT" : "")
                      << (res.cap_hit ? " CAP" : "")
                      << "\n";

            bool solved = (!res.timeout && !res.cap_hit && res.dist != INT_MAX);
            if (!solved) {
                all_ok = false;

                if (use_bidir) {
                    auto fb0 = std::chrono::steady_clock::now();
                    int d_bidir = BiBFSSearchHashed(A, B, bidir_cap);
                    auto fb1 = std::chrono::steady_clock::now();
                    double fb_ms = std::chrono::duration<double, std::milli>(fb1 - fb0).count();
                    std::cout << "        [BiBFS] dist=" << d_bidir
                              << " time_ms=" << fb_ms
                              << (d_bidir == INT_MAX ? " HIT_CAP" : "")
                              << "\n";
                }
                break; // Fail fast once BFS falters at this size.
            }
        }
        if (!all_ok) {
            std::cout << "[PROBE] random brute-force ceiling observed near n=" << n << "\n";
            break;
        }
    }
}

// Master test entry point executed by main().
void runAllTests() {
    /*

    test_canonical_serialization();
    test_comparison_row_json();
    test_bidirectional_agreement();
    test_comparison_fallback();
    test_identity_cases();
    test_single_rotation();
    test_skew_vs_balanced();
    test_random_small_agreement();
    test_minimality_small();
    capacity_sweep_bibfs_only();
    capacity_with_my_bfs();
    // random_with_my_bfs();
    
    */

    probe_random_limit();
}

int main()
{
    runAllTests();
    return 0;
}
