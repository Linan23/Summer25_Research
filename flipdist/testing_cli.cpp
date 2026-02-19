#include "testing_cli.h"
#include "algorithm.h"
#include "helpers.h"
#include "memoization.h"
#include "bf_bst.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Definition: Simple stopwatch for milliseconds and microseconds timing
// Parameters: none
// Returns: constructed timer
// Errors: none
class PerformanceTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
public:
    // Definition: Start or reset the timer
    // Parameters: none
    // Returns: nothing
    // Errors: none
    void start() { start_time = std::chrono::high_resolution_clock::now(); }

    // Definition: Elapsed time in microseconds since start
    // Parameters: none
    // Returns: elapsed microseconds
    // Errors: none
    long long getMicroseconds() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }

    // Definition: Elapsed time in milliseconds since start
    // Parameters: none
    // Returns: elapsed milliseconds
    // Errors: none
    double getMilliseconds() { return getMicroseconds() / 1000.0; }
};

// Definition: Deterministic tree shape generators for CLI inputs
// Parameters: none
// Returns: class with static generation helpers
// Errors: none
class TreeGenerator {
public:
    // Definition: Right-leaning chain tree
    // Parameters: n: number of nodes
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateRightChain(int n) {
        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++) {
            preorder.push_back(i);
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    // Definition: Left-leaning chain tree
    // Parameters: n: number of nodes
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateLeftChain(int n) {
        std::vector<int> preorder, inorder;
        for (int i = n; i >= 1; i--) {
            preorder.push_back(i);
        }
        for (int i = 1; i <= n; i++) {
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    // Definition: Balanced BST shape
    // Parameters: n: number of nodes
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateBalanced(int n) {
        if (n == 0) return {{}, {}};
        if (n == 1) return {{1}, {1}};
        if (n == 2) return {{2, 1}, {1, 2}};
        if (n == 3) return {{2, 1, 3}, {1, 2, 3}};
        if (n == 4) return {{3, 2, 1, 4}, {1, 2, 3, 4}};
        if (n == 5) return {{3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}};
        if (n == 6) return {{4, 2, 1, 3, 5, 6}, {1, 2, 3, 4, 5, 6}};
        if (n == 7) return {{4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7}};

        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++) inorder.push_back(i);

        std::function<void(int, int)> buildBalanced = [&](int start, int end) {
            if (start > end) return;
            int mid = (start + end) / 2;
            preorder.push_back(mid);
            buildBalanced(start, mid - 1);
            buildBalanced(mid + 1, end);
        };

        buildBalanced(1, n);
        return {preorder, inorder};
    }

    // Definition: Random BST shape from shuffled insertion order
    // Parameters: n: number of nodes; seed: shuffle seed
    // Returns: {preorder, inorder} vectors
    // Errors: returns empty vectors if n <= 0
    static std::pair<std::vector<int>, std::vector<int>> generateRandom(int n, int seed = 42) {
        if (n <= 0) return {{}, {}};

        struct Node {
            int v;
            Node *l = nullptr;
            Node *r = nullptr;
            explicit Node(int val) : v(val) {}
        };
        auto insert = [](Node *root, int val) -> Node * {
            if (!root) return new Node(val);
            Node *cur = root;
            while (true) {
                if (val < cur->v) {
                    if (cur->l) cur = cur->l;
                    else { cur->l = new Node(val); break; }
                } else {
                    if (cur->r) cur = cur->r;
                    else { cur->r = new Node(val); break; }
                }
            }
            return root;
        };
        auto preorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
            if (!node) return;
            out.push_back(node->v);
            self(node->l, out, self);
            self(node->r, out, self);
        };
        auto inorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
            if (!node) return;
            self(node->l, out, self);
            out.push_back(node->v);
            self(node->r, out, self);
        };

        std::vector<int> values;
        values.reserve(n);
        for (int i = 1; i <= n; i++) values.push_back(i);
        std::mt19937 rng(seed);
        std::shuffle(values.begin(), values.end(), rng);

        Node *root = nullptr;
        for (int v : values) root = insert(root, v);

        std::vector<int> pre, in;
        preorder(root, pre, preorder);
        inorder(root, in, inorder);
        return {pre, in};
    }
};

// Definition: Encode a tree as P:...;I:... for parity scripts
// Parameters: T: tree to encode
// Returns: string encoding of preorder and inorder
// Errors: returns empty encodings if the tree has empty traversals
static std::string treeToPreInString(const VectorRangeTreeMap &T) {
    auto join = [](const std::vector<int> &vals) {
        std::ostringstream oss;
        for (size_t i = 0; i < vals.size(); i++) {
            if (i) oss << ",";
            oss << vals[i];
        }
        return oss.str();
    };
    std::ostringstream oss;
    oss << "P:" << join(T.original_preorder) << ";I:" << join(T.original_inorder);
    return oss.str();
}

struct CliOptions {
    std::string case_type = "random";
    int n = 0;
    int seed = 0;
    int count = 1;
    int max_k = -1;
    int bfs_cap = 1;
    bool help = false;
};

// Definition: Emit CLI usage help
// Parameters: argv0: program name
// Returns: nothing
// Errors: none
static void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0 << " --case random|comb --n <int> [--seed <int>] [--count <int>]\n"
              << "       [--max-k <int>] [--bfs-cap <int>]\n";
}

// Definition: Parse CLI flags into CliOptions
// Parameters: argc/argv: CLI args; opts: output options
// Returns: true if parsing succeeds, else false
// Errors: prints to stderr on unknown args
static bool parseArgs(int argc, char **argv, CliOptions &opts) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            opts.help = true;
            return true;
        }
        if (arg == "--case" && i + 1 < argc) {
            opts.case_type = argv[++i];
            continue;
        }
        if (arg == "--n" && i + 1 < argc) {
            opts.n = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--seed" && i + 1 < argc) {
            opts.seed = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--count" && i + 1 < argc) {
            opts.count = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--max-k" && i + 1 < argc) {
            opts.max_k = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--bfs-cap" && i + 1 < argc) {
            opts.bfs_cap = std::stoi(argv[++i]);
            continue;
        }
        if (arg == "--emit-path" || arg == "--path-ascii") {
            continue; // accepted for compatibility
        }
        std::cerr << "Unknown argument: " << arg << "\n";
        return false;
    }
    return true;
}

// Definition: Execute the FlipDist CLI and emit JSON result lines
// Parameters: argc/argv: CLI args
// Returns: exit code, 0 on success and nonzero on invalid args
// Errors: prints to stderr on invalid args or unsupported cases
int runCli(int argc, char **argv) {
    CliOptions opts;
    if (!parseArgs(argc, argv, opts)) {
        printUsage(argv[0]);
        return 2;
    }
    if (opts.help || argc == 1) {
        printUsage(argv[0]);
        return opts.help ? 0 : 2;
    }
    if (opts.n <= 0) {
        std::cerr << "Missing or invalid --n\n";
        return 2;
    }
    if (opts.count <= 0) {
        std::cerr << "Invalid --count\n";
        return 2;
    }
    if (opts.case_type != "random" && opts.case_type != "comb") {
        std::cerr << "Invalid --case (use random or comb)\n";
        return 2;
    }
    if (opts.max_k <= 0) {
        opts.max_k = std::max(1, 3 * opts.n + 10);
    }

    for (int i = 0; i < opts.count; i++) {
        int case_seed = opts.case_type == "random" ? (opts.seed + i) : -1;

        std::pair<std::vector<int>, std::vector<int>> pre1_in1;
        std::pair<std::vector<int>, std::vector<int>> pre2_in2;
        if (opts.case_type == "random") {
            pre1_in1 = TreeGenerator::generateRandom(opts.n, case_seed * 2);
            pre2_in2 = TreeGenerator::generateRandom(opts.n, case_seed * 2 + 1);
        } else {
            pre1_in1 = TreeGenerator::generateLeftChain(opts.n);
            pre2_in2 = TreeGenerator::generateRightChain(opts.n);
        }

        VectorRangeTreeMap A, B;
        A.build(pre1_in1.first, pre1_in1.second);
        B.build(pre2_in2.first, pre2_in2.second);

        auto run_one = [&](const VectorRangeTreeMap &T1,
                           const VectorRangeTreeMap &T2,
                           const char *direction) {
            initProfile();
            if (g_profile.enabled) {
                g_profile.start_time = std::chrono::steady_clock::now();
                g_profile.calls_flipdist = 0;
                g_profile.calls_tdi = 0;
                g_profile.calls_tds = 0;
                g_profile.calls_partition = 0;
                g_profile.free_edge_hits = 0;
                g_profile.free_edge_misses = 0;
                g_profile.s_branch_calls = 0;
                g_profile.s_empty_branch_calls = 0;
                g_profile.independent_subsets_initial = 0;
                g_profile.independent_subsets_s = 0;
                g_profile.max_s_size = 0;
                g_profile.max_i_size = 0;
                g_profile.time_flipdist_ms = 0.0;
                g_profile.time_tdi_ms = 0.0;
                g_profile.time_tds_ms = 0.0;
                g_profile.time_partition_ms = 0.0;
            }
            PerformanceTimer timer;
            timer.start();
            resetMemo();
            int dist = FlipDistMinK(T1, T2, opts.max_k);
            double ms = timer.getMilliseconds();

            std::string status = (dist >= 0) ? "ok" : "not_found";
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "{\"case_type\":\"" << opts.case_type << "\""
                      << ",\"n\":" << opts.n
                      << ",\"seed\":" << case_seed
                      << ",\"direction\":\"" << direction << "\""
                      << ",\"distance\":" << dist
                      << ",\"time_ms\":" << ms
                      << ",\"status\":\"" << status << "\""
                      << ",\"tree_a\":\"" << treeToPreInString(T1) << "\""
                      << ",\"tree_b\":\"" << treeToPreInString(T2) << "\""
                      << ",\"max_k\":" << opts.max_k
                      << "}" << std::endl;

            if (g_profile.enabled) {
                std::cout << "[PROFILE] direction=" << direction << " dist=" << dist << "\n"
                          << "  - FlipDistTree calls=" << g_profile.calls_flipdist
                          << " time_ms=" << g_profile.time_flipdist_ms << "\n"
                          << "  - TreeDistI calls=" << g_profile.calls_tdi
                          << " time_ms=" << g_profile.time_tdi_ms << "\n"
                          << "  - TreeDistS calls=" << g_profile.calls_tds
                          << " time_ms=" << g_profile.time_tds_ms << "\n"
                          << "  - partition calls=" << g_profile.calls_partition
                          << " time_ms=" << g_profile.time_partition_ms << "\n"
                          << "  - free_edge hits=" << g_profile.free_edge_hits
                          << " misses=" << g_profile.free_edge_misses << "\n"
                          << "  - S-branch calls=" << g_profile.s_branch_calls
                          << " S-empty calls=" << g_profile.s_empty_branch_calls << "\n"
                          << "  - independent_subsets_initial=" << g_profile.independent_subsets_initial
                          << " max_I_size=" << g_profile.max_i_size << "\n"
                          << "  - max_S_size=" << g_profile.max_s_size << "\n";
            }
        };

        run_one(A, B, "a->b");
        run_one(B, A, "b->a");
    }

    return 0;
}
