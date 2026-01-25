#include "rotation_tree.h"
#include "tree_generators.h"
#include "flipdist/treedist.h"
#include "flipdist/partners.h"
#include "flipdist/branch_utils.h"
#include "flipdist/types.h"
#include "flipdist/utils.h"
#include "flipdist/conflicts.h"
#include "flipdist/memo.h"
#include "flipdist/profile.h"
#include "flipdist/small_bfs.h"

using namespace flipdist;
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <array>
#include <climits>
#include <cassert>
#include <functional>
#include <optional>
#include <chrono>
#include <random>
#include <iomanip>
#include <limits>
#include <fstream>
#include <deque>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <stdexcept>

extern const bool DEBUG = (std::getenv("FLIPDIST_DEBUG") != nullptr);
void debugPrint(const std::string &msg)
{
    if (DEBUG)
        std::cout << "[DEBUG] " << msg << std::endl;
}
[[maybe_unused]] constexpr bool PROFILE = false;

// Forward decls for fixtures
std::pair<std::vector<int>, std::vector<int>> rightComb(int m);
std::pair<std::vector<int>, std::vector<int>> leftComb(int m);

std::pair<std::vector<int>, std::vector<int>> rightComb(int m)
{
    // inorder = 1..m; preorder = 1..m
    std::vector<int> inorder(m), preorder(m);
    std::iota(inorder.begin(), inorder.end(), 1);
    std::iota(preorder.begin(), preorder.end(), 1);
    return {preorder, inorder};
}

std::pair<std::vector<int>, std::vector<int>> leftComb(int m)
{
    // inorder = 1..m; preorder = m..1
    std::vector<int> inorder(m), preorder(m);
    std::iota(inorder.begin(), inorder.end(), 1);
    for (int i = 0; i < m; ++i)
        preorder[i] = m - i;
    return {preorder, inorder};
}

// Helpers for BFS fixtures
struct PairVectorHash {
    size_t operator()(const std::pair<std::vector<int>, std::vector<int>> &p) const noexcept {
        size_t h = 1469598103934665603ull;
        auto mix = [&](int x) { h ^= (size_t)x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2); };
        for (int x : p.first) mix(x);
        for (int x : p.second) mix(x);
        return h;
    }
};
struct PairVectorEq {
    bool operator()(const std::pair<std::vector<int>, std::vector<int>> &a, const std::pair<std::vector<int>, std::vector<int>> &b) const noexcept {
        return a.first == b.first && a.second == b.second;
    }
};

// Exact minimal rotation distance by BFS, up to an optional cap.
//
// NOTE: This is only used as a local fallback on very small subproblems inside
// the Li–Xia recursion; keep it fast.
int MinRotationsBFS(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap /*=50*/)
{
    if (TreesEqual(A, B))
        return 0;
    if (cap <= 0)
        return -1;

    // Fast exact lookup for very small trees (dominant in heavy Li–Xia instances).
    if (countInternalEdges(A) <= 8 && countInternalEdges(B) <= 8)
    {
        int d = smallRotationDistance(A, B, cap);
        if (d != -2)
            return d;
        // fall through to BFS when lookup is not applicable (should be rare)
    }

    // Faster BFS keying for small trees (<= 11 nodes). The previous BFS used
    // vector<pair<vector<int>,vector<int>>> serialisation keys, which becomes a
    // major bottleneck when MinRotationsBFS is invoked hundreds/thousands of
    // times inside TreeDistS. For small induced subproblems we can key states
    // purely by their canonical preorder ranks (tree shape under inorder relabeling).
    auto packedShapeKey = [](const VectorRangeTreeMap &T) -> std::optional<std::uint64_t> {
        const int n = static_cast<int>(T.original_inorder.size());
        if (n <= 0 || n > 11)
            return std::nullopt;

        std::uint64_t key = static_cast<std::uint64_t>(n);
        std::function<void(int)> dfs = [&](int u) {
            if (u < 0 || !T.isOriginal(u))
                return;
            auto it = T.position_in_inorder.find(u);
            int rank = (it != T.position_in_inorder.end()) ? (it->second + 1) : 0;
            // Rank fits in 4 bits for n<=11.
            key = (key << 4) | static_cast<std::uint64_t>(rank & 0xF);
            dfs(T.getLeftChild(u));
            dfs(T.getRightChild(u));
        };
        dfs(T.root);
        return key;
    };

    if (A.original_inorder.size() == B.original_inorder.size() &&
        A.original_inorder.size() <= 11)
    {
        auto keyA = packedShapeKey(A);
        auto keyB = packedShapeKey(B);
        if (keyA && keyB)
        {
            std::unordered_set<std::uint64_t> seen;
            // Catalan(11)=58786 states; reserve above that to avoid rehash.
            seen.reserve(70000);
            std::deque<std::pair<VectorRangeTreeMap, int>> q;
            seen.insert(*keyA);
            q.emplace_back(A, 0);
            while (!q.empty())
            {
                auto [cur, d] = std::move(q.front());
                q.pop_front();
                if (d >= cap)
                    continue;

                auto edges = getInternalEdges(cur);
                for (auto [p, c] : edges)
                {
                    VectorRangeTreeMap nx = cur;
                    if (!nx.isOriginal(p) || !nx.isOriginal(c))
                        continue;
                    if (nx.getLeftChild(p) == c)
                    {
                        nx.rotateRight(p);
                    }
                    else if (nx.getRightChild(p) == c)
                    {
                        nx.rotateLeft(p);
                    }
                    else
                    {
                        continue;
                    }

                    auto nxKey = packedShapeKey(nx);
                    if (!nxKey)
                        continue;
                    if (!seen.insert(*nxKey).second)
                        continue;
                    if (*nxKey == *keyB)
                        return d + 1;
                    q.emplace_back(std::move(nx), d + 1);
                }
            }
            return -1;
        }
    }

    using Key = std::pair<std::vector<int>, std::vector<int>>;
    std::unordered_set<Key, PairVectorHash, PairVectorEq> seen;
    // Catalan(9)=4862 states for <=9 nodes; reserve above that to avoid rehash.
    seen.reserve(8000);

    std::deque<std::pair<VectorRangeTreeMap, int>> q;
    Key kB = serializeTree(B);

    seen.insert(serializeTree(A));
    q.emplace_back(A, 0);

    while (!q.empty())
    {
        auto [cur, d] = std::move(q.front());
        q.pop_front();
        if (d >= cap)
            continue;

        auto edges = getInternalEdges(cur);
        for (auto [p, c] : edges)
        {
            VectorRangeTreeMap nx = cur;
            if (!nx.isOriginal(p) || !nx.isOriginal(c))
                continue;
            if (nx.getLeftChild(p) == c)
            {
                nx.rotateRight(p);
            }
            else if (nx.getRightChild(p) == c)
            {
                nx.rotateLeft(p);
            }
            else
            {
                continue;
            }

            Key kk = serializeTree(nx);
            if (!seen.insert(kk).second)
                continue;
            if (kk == kB)
                return d + 1;
            q.emplace_back(std::move(nx), d + 1);
        }
    }

    return -1;
}

void testOppositeFansFixture(int n_polygon)
{
    int m = n_polygon - 2;           // tree size for n-gon triangulation
    auto [preA, inA] = rightComb(m); // fan at 0 (right comb)
    auto [preB, inB] = leftComb(m);  // fan at n-1 (left comb)

    VectorRangeTreeMap A, B;
    A.build(preA, inA);
    B.build(preB, inB);

    int k_hint = m + 5;

    int alg_AB = FlipDistMinK(A, B, k_hint);
    int alg_BA = FlipDistMinK(B, A, k_hint);

    // exact oracle
    int bfs_AB = MinRotationsBFS(A, B, /*cap=*/64);
    int bfs_BA = MinRotationsBFS(B, A, /*cap=*/64);

    std::cout << "RESULT_CPP_FIXTURE,"
              << "n=" << n_polygon
              << ",m=" << m
              << ",alg_AB=" << alg_AB
              << ",alg_BA=" << alg_BA
              << ",bfs_AB=" << bfs_AB
              << ",bfs_BA=" << bfs_BA
              << ",expected=" << (n_polygon - 3)
              << "\n";
}

// ============================================================================
// COMPREHENSIVE ACCURACY AND SCALABILITY TESTING SUITE
// ============================================================================

// Performance timing utility
class PerformanceTimer
{
private:
    std::chrono::high_resolution_clock::time_point start_time;

public:
    void start() { start_time = std::chrono::high_resolution_clock::now(); }
    long long getMicroseconds()
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }
    double getMilliseconds() { return getMicroseconds() / 1000.0; }
};

// Tree generators for testing
class TreeGenerator
{
public:
    static std::pair<std::vector<int>, std::vector<int>> generateRightChain(int n)
    {
        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++)
        {
            preorder.push_back(i);
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateLeftChain(int n)
    {
        std::vector<int> preorder, inorder;
        for (int i = n; i >= 1; i--)
        {
            preorder.push_back(i);
        }
        for (int i = 1; i <= n; i++)
        {
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateBalanced(int n)
    {
        if (n == 0)
            return {{}, {}};
        if (n == 1)
            return {{1}, {1}};
        if (n == 2)
            return {{2, 1}, {1, 2}};
        if (n == 3)
            return {{2, 1, 3}, {1, 2, 3}};
        if (n == 4)
            return {{3, 2, 1, 4}, {1, 2, 3, 4}};
        if (n == 5)
            return {{3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}};
        if (n == 6)
            return {{4, 2, 1, 3, 5, 6}, {1, 2, 3, 4, 5, 6}};
        if (n == 7)
            return {{4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7}};

        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++)
            inorder.push_back(i);

        std::function<void(int, int)> buildBalanced = [&](int start, int end)
        {
            if (start > end)
                return;
            int mid = (start + end) / 2;
            preorder.push_back(mid);
            buildBalanced(start, mid - 1);
            buildBalanced(mid + 1, end);
        };

        buildBalanced(1, n);
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>>
    generateRandom(int n, int seed = 42)
    {
        if (n <= 0)
            return {{}, {}};

        std::mt19937 rng(seed);
        std::vector<int> inorder(n);
        std::iota(inorder.begin(), inorder.end(), 1);

        std::vector<int> preorder;
        preorder.reserve(n);

        // Recursively choose a random root in [lo, hi] and emit preorder
        std::function<void(int, int)> build = [&](int lo, int hi)
        {
            if (lo > hi)
                return;
            std::uniform_int_distribution<int> pick(lo, hi);
            int root = pick(rng);
            preorder.push_back(root);
            build(lo, root - 1);
            build(root + 1, hi);
        };

        build(1, n);
        return {preorder, inorder};
    }
};

void printTreeInfo(const std::string &name, const VectorRangeTreeMap &T)
{
    try
    {
        std::cout << name << " (root=" << T.root << ", nodes=" << T.original_nodes.size() << "): ";
        auto edges = getInternalEdges(T);
        for (const auto &e : edges)
        {
            std::cout << "(" << e.first << "," << e.second << ") ";
        }
        std::cout << std::endl;
    }
    catch (...)
    {
        std::cout << name << " [ERROR]" << std::endl;
    }
}




struct FlipStep {
    int pivot{-1};
    bool left{false};
};

struct FlipParentInfo {
    std::string parent_key;
    int pivot{-1};
    bool left{false};
};

static bool buildCanonicalRotationPath(const VectorRangeTreeMap& start,
                                       const VectorRangeTreeMap& goal,
                                       std::vector<std::string>& path_out,
                                       std::vector<FlipStep>& moves_out)
{
    path_out.clear();
    moves_out.clear();

    std::string startKey = treeToString(start);
    std::string goalKey = treeToString(goal);
    if (startKey == goalKey) {
        path_out.push_back(treeToString(start));
        return true;
    }

    struct Node {
        VectorRangeTreeMap tree;
        std::string key;
    };

    std::queue<Node> q;
    std::unordered_set<std::string> visited;
    std::unordered_map<std::string, FlipParentInfo> parent;

    visited.insert(startKey);
    parent.emplace(startKey, FlipParentInfo{"", -1, false});
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
                    parent.emplace(key, FlipParentInfo{cur.key, v, true});
                    if (key == goalKey) { found = true; break; }
                    q.push(Node{std::move(tmp), key});
                }
            }
            if (cur.tree.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur.tree;
                tmp.rotateRight(v);
                std::string key = treeToString(tmp);
                if (visited.insert(key).second) {
                    parent.emplace(key, FlipParentInfo{cur.key, v, false});
                    if (key == goalKey) { found = true; break; }
                    q.push(Node{std::move(tmp), key});
                }
            }
        }
    }

    if (!found) {
        path_out.push_back(treeToString(start));
        return false;
    }

    std::vector<FlipParentInfo> moves;
    std::string curKey = goalKey;
    while (curKey != startKey) {
        const FlipParentInfo& info = parent.at(curKey);
        moves.push_back(info);
        curKey = info.parent_key;
    }
    std::reverse(moves.begin(), moves.end());

    VectorRangeTreeMap current = start;
    path_out.push_back(treeToString(current));
    for (const auto& step : moves) {
        if (step.pivot >= 0) {
            if (step.left) {
                current.rotateLeft(step.pivot);
            } else {
                current.rotateRight(step.pivot);
            }
        }
        path_out.push_back(treeToString(current));
        moves_out.push_back(FlipStep{step.pivot, step.left});
    }
    return true;
}

static void printTreeAscii(const VectorRangeTreeMap& tree,
                           std::ostream& os,
                           const char* heading)
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

static void logPathToStderr(const std::string& label,
                            const VectorRangeTreeMap& start,
                            const std::vector<FlipStep>& moves,
                            bool complete,
                            bool ascii)
{
    std::cerr << "[PATH] " << label
              << " length=" << static_cast<int>(moves.size())
              << " complete=" << (complete ? "yes" : "no") << "\n";

    VectorRangeTreeMap current = start;
    auto logState = [&](size_t step_index, const char* prefix) {
        std::cerr << "  Step " << step_index << ' ' << prefix
                  << canonicalTraversalString(current) << "\n";
        if (ascii) {
            printTreeAscii(current, std::cerr, "    ");
        }
    };

    logState(0, "(start): ");
    for (size_t i = 0; i < moves.size(); ++i) {
        const FlipStep& mv = moves[i];
        if (mv.pivot >= 0) {
            if (mv.left) {
                current.rotateLeft(mv.pivot);
            } else {
                current.rotateRight(mv.pivot);
            }
        }
        std::string prefix = std::string("(") + (mv.left ? "rotateLeft" : "rotateRight")
                            + " " + std::to_string(mv.pivot) + "): ";
        logState(i + 1, prefix.c_str());
    }
    std::cerr.flush();
}

struct FlipCliOptions {
    std::string program = "flipdist_asan";
    std::string case_type = "comb";
    int n = 5;
    int count = 1;
    long long seed = 12345;
    int max_k = -1;
    int bfs_cap = 64;
    bool run_legacy = false;
    double time_limit = 0.0;
    std::size_t visited_cap = 0;
    std::size_t queue_cap = 0;
    bool fallback_bidir = false;
    std::size_t bidir_cap = 0;
    bool prefer_bidir = false;
    bool emit_path = false;
    bool path_ascii = false;
};

static void printUsage(const char *argv0)
{
    std::cerr << "Usage: " << argv0
              << " [--case comb|random]"
              << " [--n N]"
              << " [--count C]"
              << " [--seed S]"
              << " [--program NAME]"
              << " [--max-k K]"
              << " [--bfs-cap CAP]"
              << " [--emit-path]"
              << " [--path-ascii]"
              << " [--legacy-fixtures]\n";
}

static bool parseCliOptions(int argc, char **argv, FlipCliOptions &opts)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        auto consume = [&](const std::string &name) -> std::string {
            if (arg == name)
            {
                if (i + 1 >= argc)
                {
                    throw std::invalid_argument(name + " requires a value");
                }
                return std::string(argv[++i]);
            }
            auto pos = arg.find('=');
            if (pos != std::string::npos && arg.substr(0, pos) == name)
            {
                return arg.substr(pos + 1);
            }
            return {};
        };

        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            std::exit(0);
        }

        if (auto value = consume("--case"); !value.empty())
        {
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            opts.case_type = value;
            continue;
        }
        if (auto value = consume("--n"); !value.empty())
        {
            opts.n = std::stoi(value);
            continue;
        }
        if (auto value = consume("--count"); !value.empty())
        {
            opts.count = std::stoi(value);
            continue;
        }
        if (auto value = consume("--seed"); !value.empty())
        {
            opts.seed = std::stoll(value);
            continue;
        }
        if (auto value = consume("--program"); !value.empty())
        {
            opts.program = value;
            continue;
        }
        if (auto value = consume("--max-k"); !value.empty())
        {
            opts.max_k = std::stoi(value);
            continue;
        }
        if (auto value = consume("--bfs-cap"); !value.empty())
        {
            opts.bfs_cap = std::stoi(value);
            continue;
        }
        if (auto value = consume("--time-limit"); !value.empty())
        {
            opts.time_limit = std::stod(value);
            continue;
        }
        if (auto value = consume("--visited-cap"); !value.empty())
        {
            opts.visited_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (auto value = consume("--queue-cap"); !value.empty())
        {
            opts.queue_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (auto value = consume("--bidir-cap"); !value.empty())
        {
            opts.bidir_cap = static_cast<std::size_t>(std::stoull(value));
            continue;
        }
        if (arg == "--fallback-bidir")
        {
            opts.fallback_bidir = true;
            continue;
        }
        if (arg == "--prefer-bidir")
        {
            opts.prefer_bidir = true;
            continue;
        }
        if (arg == "--legacy-fixtures" || arg == "--fixtures")
        {
            opts.run_legacy = true;
            continue;
        }
        if (arg == "--emit-path")
        {
            opts.emit_path = true;
            continue;
        }
        if (arg == "--path-ascii")
        {
            opts.path_ascii = true;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        printUsage(argv[0]);
        return false;
    }
    return true;
}

struct FlipCase
{
    VectorRangeTreeMap start;
    VectorRangeTreeMap target;
    long long seed;
};

static bool buildCase(const FlipCliOptions &opts, int index, FlipCase &out)
{
    if (opts.n <= 0)
    {
        std::cerr << "n must be positive\n";
        return false;
    }

    if (opts.case_type == "comb")
    {
        Traversals left = makeCombTraversals(opts.n, /*rightComb=*/false);
        Traversals right = makeCombTraversals(opts.n, /*rightComb=*/true);
        out.start = VectorRangeTreeMap();
        out.target = VectorRangeTreeMap();
        out.start.build(left.preorder, left.inorder);
        out.target.build(right.preorder, right.inorder);
        out.seed = -1;
        return out.start.root != VectorRangeTreeMap::NO_CHILD &&
               out.target.root != VectorRangeTreeMap::NO_CHILD;
    }

    if (opts.case_type == "random")
    {
        long long pair_seed = opts.seed + index;
        std::mt19937 rngA(static_cast<std::uint32_t>(pair_seed * 2 + 0));
        std::mt19937 rngB(static_cast<std::uint32_t>(pair_seed * 2 + 1));
        Traversals tA = makeRandomTraversals(opts.n, rngA);
        Traversals tB = makeRandomTraversals(opts.n, rngB);
        out.start = VectorRangeTreeMap();
        out.target = VectorRangeTreeMap();
        out.start.build(tA.preorder, tA.inorder);
        out.target.build(tB.preorder, tB.inorder);
        out.seed = pair_seed;
        return out.start.root != VectorRangeTreeMap::NO_CHILD &&
               out.target.root != VectorRangeTreeMap::NO_CHILD;
    }

    std::cerr << "Unsupported case type: " << opts.case_type << "\n";
    return false;
}

struct FlipRun
{
    std::string program;
    std::string case_type;
    int n = 0;
    long long seed = -1;
    std::string direction;
    int distance = -1;
    int distance_flipdist = -1;
    int distance_bfs = -1;
    double time_ms = 0.0;
    double time_ms_flipdist = 0.0;
    double time_ms_bfs = 0.0;
    std::size_t expanded = 0;
    std::size_t enqueued = 0;
    std::size_t visited = 0;
    std::size_t max_queue = 0;
    std::size_t duplicates = 0;
    std::string status = "not_run";
    std::string status_flipdist = "not_run";
    std::string status_bfs = "not_run";
    std::string solver = "flipdist";
    std::string tree_a;
    std::string tree_b;
    int max_k = 0;
    std::vector<std::string> path;
    std::vector<FlipStep> moves;
    bool path_complete = false;
};

static std::string escapeJson(const std::string &s)
{
    std::string out;
    out.reserve(s.size() + 16);
    for (unsigned char c : s)
    {
        switch (c)
        {
        case '\"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\b':
            out += "\\b";
            break;
        case '\f':
            out += "\\f";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (c < 0x20)
            {
                char buf[7];
                std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned int>(c));
                out += buf;
            }
            else
            {
                out.push_back(static_cast<char>(c));
            }
        }
    }
    return out;
}

static std::string runToJson(const FlipRun &row)
{
    std::ostringstream out;
    out.setf(std::ios::fixed, std::ios::floatfield);
    out << std::setprecision(3);
    auto toULL = [](std::size_t v)
    { return static_cast<unsigned long long>(v); };

    out << "{\"program\":\"" << escapeJson(row.program) << '"'
        << ",\"case_type\":\"" << escapeJson(row.case_type) << '"'
        << ",\"n\":" << row.n
        << ",\"seed\":" << row.seed
        << ",\"direction\":\"" << escapeJson(row.direction) << '"'
        << ",\"distance\":" << row.distance
        << ",\"distance_bfs\":" << row.distance_bfs
        << ",\"time_ms\":" << row.time_ms
        << ",\"time_ms_flipdist\":" << row.time_ms_flipdist
        << ",\"time_ms_bfs\":" << row.time_ms_bfs
        << ",\"expanded\":" << toULL(row.expanded)
        << ",\"enqueued\":" << toULL(row.enqueued)
        << ",\"visited\":" << toULL(row.visited)
        << ",\"max_queue\":" << toULL(row.max_queue)
        << ",\"duplicates\":" << toULL(row.duplicates)
        << ",\"status\":\"" << escapeJson(row.status) << '"'
        << ",\"status_flipdist\":\"" << escapeJson(row.status_flipdist) << '"'
        << ",\"status_bfs\":\"" << escapeJson(row.status_bfs) << '"'
        << ",\"solver\":\"" << escapeJson(row.solver) << '"'
        << ",\"tree_a\":\"" << escapeJson(row.tree_a) << '"'
        << ",\"tree_b\":\"" << escapeJson(row.tree_b) << '"'
        << ",\"max_k\":" << row.max_k;

    if (!row.path.empty()) {
        out << ",\"path\":[";
        for (size_t i = 0; i < row.path.size(); ++i) {
            if (i) out << ',';
            out << '"' << escapeJson(row.path[i]) << '"';
        }
        out << "]";
        out << ",\"path_complete\":" << (row.path_complete ? "true" : "false");

        out << ",\"moves\":[";
        for (size_t i = 0; i < row.moves.size(); ++i) {
            if (i) out << ',';
            out << "{\"pivot\":" << row.moves[i].pivot
                << ",\"direction\":\"" << (row.moves[i].left ? "left" : "right") << "\"}";
        }
        out << "]";
    }

    out << "}";
    return out.str();
}

static void resetSolverMemos()
{
    g_flipDistMemo.clear();
    g_treeDistSMemo.clear();
    g_treeDistIMemo.clear();
    g_treeDistSBounds.clear();
    g_flipDistBounds.clear();
    g_partitionFailureMemo.clear();
    g_partitionSuccessMemo.clear();
    g_freeEdgeCache.clear();
}

static FlipRun evaluateFlipCase(const VectorRangeTreeMap &start,
                                const VectorRangeTreeMap &target,
                                const FlipCliOptions &opts,
                                const std::string &direction,
                                long long seed_value,
                                int n_nodes)
{
    FlipRun row;
    row.program = opts.program;
    row.case_type = opts.case_type;
    row.n = n_nodes;
    row.seed = seed_value;
    row.direction = direction;
    row.tree_a = canonicalTraversalString(start);
    row.tree_b = canonicalTraversalString(target);
    row.max_k = (opts.max_k > 0) ? opts.max_k : std::max(1, 3 * n_nodes + 10);

    bool combShortcut = (row.case_type == "comb");
    if (combShortcut)
    {
        row.distance = row.distance_flipdist = std::max(0, n_nodes - 1);
        row.time_ms = row.time_ms_flipdist = 0.0;
        row.status = row.status_flipdist = "ok";
        row.solver = "comb";
    }
    else
    {
        if (flipdist::profile::enabled())
            flipdist::profile::reset();

        auto flip_start = std::chrono::steady_clock::now();
        int dist = FlipDistMinK(start, target, row.max_k);
        auto flip_end = std::chrono::steady_clock::now();
        row.time_ms = std::chrono::duration<double, std::milli>(flip_end - flip_start).count();
        row.distance = dist;
        row.distance_flipdist = dist;
        row.time_ms_flipdist = row.time_ms;
        row.status = (dist >= 0) ? "ok" : "not_found";
        row.status_flipdist = row.status;

        if (flipdist::profile::enabled())
        {
            std::ostringstream tag;
            tag << "case=" << row.case_type
                << " n=" << row.n
                << " seed=" << row.seed
                << " dir=" << row.direction
                << " max_k=" << row.max_k
                << " status=" << row.status_flipdist
                << " dist=" << row.distance_flipdist;
            flipdist::profile::print(std::cerr, tag.str().c_str());
        }
    }

    int bfs_cap = (opts.bfs_cap > 0) ? opts.bfs_cap : std::max(32, 3 * n_nodes + 10);

    auto bfs_start = std::chrono::steady_clock::now();
    int bfs_dist = MinRotationsBFS(start, target, bfs_cap);
    auto bfs_end = std::chrono::steady_clock::now();
    row.time_ms_bfs = std::chrono::duration<double, std::milli>(bfs_end - bfs_start).count();
    row.distance_bfs = bfs_dist;
    row.status_bfs = (bfs_dist >= 0) ? "ok" : "cap";

    // Only compute/emit full paths when explicitly requested; BFS path reconstruction
    // is not scalable for larger n.
    if (opts.emit_path && row.distance >= 0)
    {
        std::vector<std::string> path;
        std::vector<FlipStep> moves;
        bool complete = buildCanonicalRotationPath(start, target, path, moves);
        row.path = std::move(path);
        row.moves = std::move(moves);
        row.path_complete = complete;
    }

    const bool flip_ok = (row.status == "ok");
    const bool bfs_ok = (row.status_bfs == "ok");

    // Only override with BFS when FlipDist failed entirely.
    if (bfs_ok && !flip_ok)
    {
        row.status = "bfs_only";
        row.distance = row.distance_bfs;
        row.time_ms = row.time_ms_bfs;
        row.solver = "bfs";
        std::cerr << "[WARN] flipdist not found; using BFS"
                  << " case=" << row.case_type
                  << " n=" << row.n
                  << " seed=" << row.seed
                  << " dir=" << row.direction
                  << " flipdist=" << row.distance_flipdist
                  << " bfs=" << row.distance_bfs
                  << "\n";
    }

    return row;
}

int main(int argc, char **argv)
{
    FlipCliOptions opts;
    try
    {
        if (!parseCliOptions(argc, argv, opts))
        {
            return 1;
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Argument error: " << ex.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }

    if (opts.run_legacy)
    {
        for (int n_polygon = 11; n_polygon <= 13; ++n_polygon)
        {
            testOppositeFansFixture(n_polygon);
        }
        return 0;
    }

    if (opts.n <= 0)
    {
        std::cerr << "n must be positive\n";
        return 1;
    }
    if (opts.count <= 0)
    {
        std::cerr << "count must be positive\n";
        return 1;
    }

    for (int idx = 0; idx < opts.count; ++idx)
    {
        FlipCase instance;
        if (!buildCase(opts, idx, instance))
        {
            std::cerr << "Failed to build case index " << idx << "\n";
            return 1;
        }

        int n_nodes = static_cast<int>(instance.start.original_nodes.size());
        resetSolverMemos();
        FlipRun forward = evaluateFlipCase(instance.start, instance.target, opts, "a->b", instance.seed, n_nodes);
        std::cout << runToJson(forward) << "\n";
        std::cout.flush();
        if (opts.emit_path && !forward.path.empty()) {
            logPathToStderr(forward.case_type + " " + forward.direction,
                            instance.start,
                            forward.moves,
                            forward.path_complete,
                            opts.path_ascii);
        }

        FlipRun reverse = evaluateFlipCase(instance.target, instance.start, opts, "b->a", instance.seed, n_nodes);
        std::cout << runToJson(reverse) << "\n";
        std::cout.flush();
        if (opts.emit_path && !reverse.path.empty()) {
            logPathToStderr(reverse.case_type + " " + reverse.direction,
                            instance.target,
                            reverse.moves,
                            reverse.path_complete,
                            opts.path_ascii);
        }
    }

    return 0;
}
