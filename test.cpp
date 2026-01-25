// Houses the lightweight CLI harness we use during development. It keeps the
// legacy unit checks, the random/comb scaling probes, and the interactive demo
// (`BFS_DEMO=1`) that prints trees and reconstructs the exact rotation path
// in both directions. All behaviour is controlled through environment
// variables so we can trigger probes, cap sizes, or request visualisations
// without recompiling.

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
#include <functional>

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

// Verifies that partitionAlongEdge reconstructs the actual current subtree
// structure (not just the node set) after arbitrary rotations.
static void test_partition_along_edge_consistency()
{
    if (!envAsFlag("TEST_PARTITION_ALONG_EDGE"))
        return;

    const int n = static_cast<int>(envAsSize("TEST_PARTITION_N", 15));
    const int rotations = static_cast<int>(envAsSize("TEST_PARTITION_ROTATIONS", 200));
    const int trials = static_cast<int>(envAsSize("TEST_PARTITION_TRIALS", 200));
    const uint32_t seed = static_cast<uint32_t>(envAsSize("TEST_PARTITION_SEED", 1));

    std::mt19937 rng(seed);
    Traversals base = makeRandomTraversals(n, rng);
    VectorRangeTreeMap T;
    T.build(base.preorder, base.inorder);

    auto randomRotateOnce = [&](VectorRangeTreeMap &tree) {
        std::vector<int> pivots;
        pivots.reserve(tree.original_nodes.size());
        for (int v : tree.original_nodes)
        {
            if (tree.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD ||
                tree.getRightChild(v) != VectorRangeTreeMap::NO_CHILD)
                pivots.push_back(v);
        }
        if (pivots.empty())
            return;

        std::uniform_int_distribution<size_t> pickPivot(0, pivots.size() - 1);
        int pivot = pivots[pickPivot(rng)];

        bool canLeft = (tree.getRightChild(pivot) != VectorRangeTreeMap::NO_CHILD);
        bool canRight = (tree.getLeftChild(pivot) != VectorRangeTreeMap::NO_CHILD);
        if (!canLeft && !canRight)
            return;

        bool doLeft = false;
        if (canLeft && canRight)
        {
            std::uniform_int_distribution<int> pickDir(0, 1);
            doLeft = (pickDir(rng) == 1);
        }
        else
        {
            doLeft = canLeft;
        }

        if (doLeft)
            tree.rotateLeft(pivot);
        else
            tree.rotateRight(pivot);
    };

    for (int i = 0; i < rotations; ++i)
        randomRotateOnce(T);

    auto collectDirectedEdges = [](const VectorRangeTreeMap &tree) {
        std::unordered_set<std::pair<int,int>, PairHash, PairEq> edges;
        edges.reserve(tree.original_nodes.size() * 2);
        tree.collectEdges(tree.root, edges);
        std::vector<std::pair<int,int>> out(edges.begin(), edges.end());
        std::sort(out.begin(), out.end());
        return out;
    };

    auto buildSubtreeFromNode = [](const VectorRangeTreeMap &tree, int rootNode) {
        std::vector<int> preorder;
        std::vector<int> inorder;
        preorder.reserve(tree.original_nodes.size());
        inorder.reserve(tree.original_nodes.size());

        std::function<void(int)> dfsPre = [&](int node) {
            if (node < 0 || !tree.isOriginal(node))
                return;
            preorder.push_back(node);
            dfsPre(tree.getLeftChild(node));
            dfsPre(tree.getRightChild(node));
        };
        std::function<void(int)> dfsIn = [&](int node) {
            if (node < 0 || !tree.isOriginal(node))
                return;
            dfsIn(tree.getLeftChild(node));
            inorder.push_back(node);
            dfsIn(tree.getRightChild(node));
        };

        dfsPre(rootNode);
        dfsIn(rootNode);

        VectorRangeTreeMap out;
        out.build(preorder, inorder);
        return out;
    };

    for (int t = 0; t < trials; ++t)
    {
        auto edges = collectDirectedEdges(T);
        if (edges.empty())
            break;
        std::uniform_int_distribution<size_t> pickEdge(0, edges.size() - 1);
        auto [parent, child] = edges[pickEdge(rng)];

        auto parent_range = T.getRange(parent);
        auto child_range = T.getRange(child);
        auto [A, B] = VectorRangeTreeMap::partitionAlongEdge(T, parent_range, child_range);

        VectorRangeTreeMap expected = buildSubtreeFromNode(T, child);

        if (!TreesEqual(A, expected))
        {
            std::cerr << "[FAIL] partitionAlongEdge mismatch on trial " << t
                      << " edge=(" << parent << "," << child << ")\n";
            std::cerr << "  expected=" << treeToString(expected) << "\n";
            std::cerr << "  got     =" << treeToString(A) << "\n";
            std::abort();
        }
    }

    std::cout << "[OK] partitionAlongEdge structure matches current subtree (" << trials << " trials)\n";
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

static std::string envAsString(const char* name, const std::string& fallback)
{
    if (const char* val = std::getenv(name)) return std::string(val);
    return fallback;
}

static void printTreeAscii(const VectorRangeTreeMap& tree, std::ostream& os, const char* heading)
{
    os << heading << "\n";
    if (tree.root == VectorRangeTreeMap::NO_CHILD) {
        os << "  (empty)\n";
        return;
    }

    std::function<void(int,const std::string&,bool)> rec =
        [&](int node, const std::string& prefix, bool left) {
            if (node == VectorRangeTreeMap::NO_CHILD || !tree.isOriginal(node)) return;

            int right = tree.getRightChild(node);
            int leftChild = tree.getLeftChild(node);
            bool hasRight = (right != VectorRangeTreeMap::NO_CHILD) && tree.isOriginal(right);
            bool hasLeft  = (leftChild != VectorRangeTreeMap::NO_CHILD) && tree.isOriginal(leftChild);

            if (hasRight) rec(right, prefix + (left ? "│   " : "    "), false);

            os << prefix;
            if (!prefix.empty()) {
                os << (left ? "└── " : "┌── ");
            }
            os << node << "\n";

            if (hasLeft)  rec(leftChild, prefix + (left ? "    " : "│   "), true);
        };

    rec(tree.root, "", true);
}

struct DemoParentInfo {
    std::string parent_key;
    int pivot{-1};
    bool left{false};
};

static std::vector<VectorRangeTreeMap>
buildDemoPath(const VectorRangeTreeMap& start,
              const VectorRangeTreeMap& goal,
              std::vector<DemoParentInfo>& moves_out)
{
    moves_out.clear();
    std::string startKey = treeToString(start);
    std::string goalKey  = treeToString(goal);
    if (startKey == goalKey) return {start};

    struct Node { VectorRangeTreeMap tree; std::string key; };

    std::queue<Node> q;
    std::unordered_set<std::string> visited;
    std::unordered_map<std::string, DemoParentInfo> parent;

    visited.insert(startKey);
    parent.emplace(startKey, DemoParentInfo{"", -1, false});
    q.push(Node{start, startKey});

    bool found = false;
    while (!q.empty() && !found) {
        Node cur = std::move(q.front());
        q.pop();

        for (int v : cur.tree.original_nodes) {
            if (cur.tree.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur.tree;
                tmp.rotateLeft(v);
                std::string key = treeToString(tmp);
                if (visited.insert(key).second) {
                    parent.emplace(key, DemoParentInfo{cur.key, v, true});
                    if (key == goalKey) { found = true; break; }
                    q.push(Node{tmp, key});
                }
            }
            if (cur.tree.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur.tree;
                tmp.rotateRight(v);
                std::string key = treeToString(tmp);
                if (visited.insert(key).second) {
                    parent.emplace(key, DemoParentInfo{cur.key, v, false});
                    if (key == goalKey) { found = true; break; }
                    q.push(Node{tmp, key});
                }
            }
            if (found) break;
        }
    }

    if (!found) return {start};

    std::vector<DemoParentInfo> moves;
    std::string curKey = goalKey;
    while (curKey != startKey) {
        const DemoParentInfo& info = parent.at(curKey);
        moves.push_back(info);
        curKey = info.parent_key;
    }
    std::reverse(moves.begin(), moves.end());
    moves_out = moves;

    std::vector<VectorRangeTreeMap> path;
    path.reserve(moves.size() + 1);
    VectorRangeTreeMap current = start;
    path.push_back(current);
    for (const auto& step : moves) {
        if (step.pivot >= 0) {
            if (step.left) current.rotateLeft(step.pivot);
            else           current.rotateRight(step.pivot);
        }
        path.push_back(current);
    }
    return path;
}

static void printDemoPath(const char* label,
                          const std::vector<VectorRangeTreeMap>& path,
                          const std::vector<DemoParentInfo>& moves,
                          bool ascii,
                          size_t max_steps)
{
    if (path.empty()) return;
    std::cout << label << " path length=" << (path.size() - 1) << "\n";
    size_t limit = max_steps ? std::min(max_steps, path.size()) : path.size();
    for (size_t i = 0; i < limit; ++i) {
        if (i == 0) {
            std::cout << "  Step 0 (start): ";
        } else {
            const DemoParentInfo& mv = moves[i-1];
            std::cout << "  Step " << i << " (" << (mv.left ? "rotateLeft" : "rotateRight")
                      << ' ' << mv.pivot << "): ";
        }
        std::cout << canonicalTraversalString(path[i]) << "\n";
        if (ascii) {
            printTreeAscii(path[i], std::cout, "    ");
        }
    }
    if (limit < path.size()) {
        std::cout << "  ... (" << (path.size() - limit) << " more steps omitted)\n";
    }
}

static void run_bruteforce_demo()
{
    if (!envAsFlag("BFS_DEMO")) return;

    std::string mode = envAsString("BFS_DEMO_MODE", "random");
    int n = static_cast<int>(envAsSize("BFS_DEMO_N", 10));
    uint32_t seedA = static_cast<uint32_t>(envAsSize("BFS_DEMO_SEED_A", 2025));
    uint32_t seedB = static_cast<uint32_t>(envAsSize("BFS_DEMO_SEED_B", seedA + 1));
    bool show_ascii = envAsFlag("BFS_DEMO_ASCII", true);
    size_t max_steps = envAsSize("BFS_DEMO_MAX_STEPS", 0);

    VectorRangeTreeMap A, B;
    if (mode == "comb") {
        Traversals left  = makeCombTraversals(n, /*rightComb=*/false);
        Traversals right = makeCombTraversals(n, /*rightComb=*/true);
        A.build(left.preorder, left.inorder);
        B.build(right.preorder, right.inorder);
    } else {
        std::mt19937 rngA(seedA);
        std::mt19937 rngB(seedB);
        Traversals tA = makeRandomTraversals(n, rngA);
        Traversals tB = makeRandomTraversals(n, rngB);
        A.build(tA.preorder, tA.inorder);
        B.build(tB.preorder, tB.inorder);
    }

    std::cout << "=== C++ brute-force demo ===\n";
    std::cout << "mode=" << mode << " n=" << n
              << " seedA=" << seedA << " seedB=" << seedB << "\n";
    std::cout << "Tree A canonical: " << canonicalTraversalString(A) << "\n";
    if (show_ascii) printTreeAscii(A, std::cout, "Tree A");
    std::cout << "Tree B canonical: " << canonicalTraversalString(B) << "\n";
    if (show_ascii) printTreeAscii(B, std::cout, "Tree B");

    BFSStats statsAB{}, statsBA{};
    int distAB = BFSSearchOptimized(A, B, &statsAB);
    int distBA = BFSSearchOptimized(B, A, &statsBA);

    auto dumpStats = [](const char* label, int dist, const BFSStats& st) {
        std::cout << label << " dist=" << dist
                  << " generated=" << st.generated
                  << " enqueued=" << st.enqueued
                  << " duplicates=" << st.duplicates
                  << " visited=" << st.visited
                  << " max_queue=" << st.max_queue
                  << " equality_hits=" << st.equality_hits
                  << " signature_mismatches=" << st.signature_mismatches
                  << "\n";
    };

    dumpStats("A -> B", distAB, statsAB);
    dumpStats("B -> A", distBA, statsBA);

    if (distAB == distBA) {
        std::cout << "[OK] symmetric distances match." << "\n";
    } else {
        std::cout << "[WARN] distances differ (A->B=" << distAB
                  << ", B->A=" << distBA << ")" << "\n";
    }
    std::cout.flush();

    std::vector<DemoParentInfo> movesAB, movesBA;
    auto pathAB = buildDemoPath(A, B, movesAB);
    auto pathBA = buildDemoPath(B, A, movesBA);
    printDemoPath("A -> B", pathAB, movesAB, show_ascii, max_steps);
    printDemoPath("B -> A", pathBA, movesBA, show_ascii, max_steps);
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

        const double time_limit = envAsDouble("BFS_TIME_LIMIT", 3.0);
        const size_t visited_cap = envAsSize("BFS_VISITED_CAP", 5'000'000);
        const size_t queue_cap   = envAsSize("BFS_QUEUE_CAP",   5'000'000);
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

            const double time_limit = envAsDouble("BFS_RANDOM_TIME_LIMIT", 2.0);
            const size_t visited_cap = envAsSize("BFS_RANDOM_VISITED_CAP", 5'000'000);
            const size_t queue_cap   = envAsSize("BFS_RANDOM_QUEUE_CAP",   5'000'000);
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

            if (envAsFlag("BFS_RANDOM_USE_BIDIR"))
            {
                size_t state_cap = envAsSize("BFS_RANDOM_BIDIR_CAP", 5'000'000);
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
    if (!envAsFlag("BFS_PROBE_RANDOM_LIMIT"))
        return;

    int start_n = static_cast<int>(envAsSize("BFS_PROBE_RANDOM_START", 8));
    int max_n   = static_cast<int>(envAsSize("BFS_PROBE_RANDOM_MAX",   30));
    int trials  = static_cast<int>(envAsSize("BFS_PROBE_RANDOM_TRIALS", 5));
    uint32_t seed = static_cast<uint32_t>(envAsSize("BFS_PROBE_RANDOM_SEED", 2025));

    double time_limit = envAsDouble("BFS_PROBE_RANDOM_TIME_LIMIT", 5.0);
    size_t visited_cap = envAsSize("BFS_PROBE_RANDOM_VISITED_CAP", 5'000'000);
    size_t queue_cap   = envAsSize("BFS_PROBE_RANDOM_QUEUE_CAP",   5'000'000);

    bool use_bidir = envAsFlag("BFS_PROBE_RANDOM_USE_BIDIR");
    size_t bidir_cap = envAsSize("BFS_PROBE_RANDOM_BIDIR_CAP", 5'000'000);

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

    probe_random_limit();
    run_bruteforce_demo();
    test_partition_along_edge_consistency();
}

int main()
{
    // Flush std::cout on every write so test output is visible even when piped.
    std::cout.setf(std::ios::unitbuf);
    runAllTests();
    return 0;
}
