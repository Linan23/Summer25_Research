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
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

// Definition: Encode current tree state as P:...;I:... by traversing live structure
// Parameters: T: tree to encode
// Returns: string encoding of current preorder and inorder traversals
// Errors: returns empty traversals for invalid/empty trees
static std::string currentTreeToPreInString(const VectorRangeTreeMap &T) {
    std::vector<int> preorder;
    std::vector<int> inorder;
    std::function<void(int)> dfs_pre = [&](int node) {
        if (node < 0 || !T.isOriginal(node)) return;
        preorder.push_back(node);
        dfs_pre(T.getLeftChild(node));
        dfs_pre(T.getRightChild(node));
    };
    std::function<void(int)> dfs_in = [&](int node) {
        if (node < 0 || !T.isOriginal(node)) return;
        dfs_in(T.getLeftChild(node));
        inorder.push_back(node);
        dfs_in(T.getRightChild(node));
    };
    if (T.root >= 0 && T.isOriginal(T.root)) {
        dfs_pre(T.root);
        dfs_in(T.root);
    }
    auto join = [](const std::vector<int> &vals) {
        std::ostringstream oss;
        for (size_t i = 0; i < vals.size(); i++) {
            if (i) oss << ",";
            oss << vals[i];
        }
        return oss.str();
    };
    std::ostringstream oss;
    oss << "P:" << join(preorder) << ";I:" << join(inorder);
    return oss.str();
}

// Definition: Render a horizontal ASCII tree (parent over children) similar to textbook diagrams
// Parameters: T: tree to render; base_indent: prefix added to each line
// Returns: nothing
// Errors: prints a placeholder for empty trees
static void printAsciiTreeDiagram(const VectorRangeTreeMap &T, const std::string &base_indent = "    ") {
    if (T.root < 0 || !T.isOriginal(T.root) || T.original_nodes.empty()) {
        std::cout << base_indent << "(empty tree)\n";
        return;
    }

    std::vector<int> inorder_nodes;
    std::unordered_map<int, int> depth_by_node;
    int max_depth = 0;

    std::function<void(int)> dfs_in = [&](int node) {
        if (node < 0 || !T.isOriginal(node)) return;
        dfs_in(T.getLeftChild(node));
        inorder_nodes.push_back(node);
        dfs_in(T.getRightChild(node));
    };
    std::function<void(int, int)> dfs_depth = [&](int node, int depth) {
        if (node < 0 || !T.isOriginal(node)) return;
        depth_by_node[node] = depth;
        if (depth > max_depth) max_depth = depth;
        dfs_depth(T.getLeftChild(node), depth + 1);
        dfs_depth(T.getRightChild(node), depth + 1);
    };

    dfs_in(T.root);
    dfs_depth(T.root, 0);

    std::unordered_map<int, int> rank_by_node;
    for (size_t i = 0; i < inorder_nodes.size(); i++) {
        rank_by_node[inorder_nodes[i]] = static_cast<int>(i);
    }

    const int spacing = 5;
    const int width = std::max(8, static_cast<int>(inorder_nodes.size()) * spacing + 2);
    const int height = max_depth * 2 + 1;
    std::vector<std::string> canvas(height, std::string(width, ' '));

    auto safe_set = [&](int r, int c, char ch) {
        if (r < 0 || r >= height || c < 0 || c >= width) return;
        canvas[r][c] = ch;
    };
    auto draw_hline = [&](int r, int c1, int c2, char ch) {
        if (r < 0 || r >= height) return;
        if (c1 > c2) std::swap(c1, c2);
        c1 = std::max(c1, 0);
        c2 = std::min(c2, width - 1);
        for (int c = c1; c <= c2; c++) canvas[r][c] = ch;
    };
    auto center_x = [&](int node) {
        return rank_by_node[node] * spacing + spacing / 2;
    };
    auto put_label = [&](int r, int x, const std::string &label) {
        int start = x - static_cast<int>(label.size()) / 2;
        for (size_t i = 0; i < label.size(); i++) {
            int c = start + static_cast<int>(i);
            if (r >= 0 && r < height && c >= 0 && c < width) {
                canvas[r][c] = label[i];
            }
        }
    };

    for (int node : inorder_nodes) {
        int y = depth_by_node[node] * 2;
        int x = center_x(node);
        put_label(y, x, std::to_string(node));

        int l = T.getLeftChild(node);
        if (l >= 0 && T.isOriginal(l)) {
            int cx = center_x(l);
            int row = y + 1;
            safe_set(row, cx + 1, '/');
            draw_hline(row, cx + 2, x - 1, '_');
        }
        int r = T.getRightChild(node);
        if (r >= 0 && T.isOriginal(r)) {
            int cx = center_x(r);
            int row = y + 1;
            safe_set(row, cx - 1, '\\');
            draw_hline(row, x + 1, cx - 2, '_');
        }
    }

    for (const auto &raw : canvas) {
        size_t end = raw.find_last_not_of(' ');
        if (end == std::string::npos) continue;
        std::cout << base_indent << raw.substr(0, end + 1) << "\n";
    }
}

// Definition: Print node relationships for a tree state
// Parameters: T: tree to print; indent: line prefix
// Returns: nothing
// Errors: prints nothing for empty trees
static void printNodeRelations(const VectorRangeTreeMap &T, const std::string &indent = "    ") {
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());
    if (nodes.empty()) return;

    std::cout << indent << "node parent left right range\n";
    for (int v : nodes) {
        if (!T.isOriginal(v)) continue;
        int p = T.getParent(v);
        int l = T.getLeftChild(v);
        int r = T.getRightChild(v);
        auto rg = T.getRange(v);
        std::cout << indent << v << " "
                  << p << " "
                  << l << " "
                  << r << " ["
                  << rg.first << "," << rg.second << "]\n";
    }
}

struct CliOptions {
    std::string case_type = "random";
    int n = 0;
    int seed = 0;
    int count = 1;
    int max_k = -1;
    int bfs_cap = 1;
    bool print_trees = false;
    bool emit_path = false;
    bool path_ascii = false;
    bool help = false;
};

// Definition: Emit CLI usage help
// Parameters: argv0: program name
// Returns: nothing
// Errors: none
static void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0 << " --case random|comb --n <int> [--seed <int>] [--count <int>]\n"
              << "       [--max-k <int>] [--bfs-cap <int>] [--print-trees] [--emit-path] [--path-ascii]\n";
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
        if (arg == "--print-trees") {
            opts.print_trees = true;
            continue;
        }
        if (arg == "--emit-path") {
            opts.emit_path = true;
            opts.path_ascii = true;
            continue;
        }
        if (arg == "--path-ascii") {
            opts.path_ascii = true;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << "\n";
        return false;
    }
    return true;
}

struct PathParent {
    std::string parent_key;
    std::string move;
    int depth = 0;
};

struct PathBuildResult {
    bool found = false;
    bool capped = false;
    std::vector<VectorRangeTreeMap> states;
    std::vector<std::string> moves;
};

// Definition: Build a shortest rotation path with BFS for debug printing
// Parameters: source/target: trees; depth_cap: max BFS depth; state_cap: max visited states
// Returns: path states and moves if found, or capped/not-found status
// Errors: returns found=false on cap hit or if no path within depth_cap
static PathBuildResult buildShortestRotationPathBfs(const VectorRangeTreeMap &source,
                                                    const VectorRangeTreeMap &target,
                                                    int depth_cap,
                                                    size_t state_cap) {
    PathBuildResult out;
    if (depth_cap < 0) return out;

    const std::string src_key = treeToString(source);
    const std::string dst_key = treeToString(target);
    if (src_key == dst_key) {
        out.found = true;
        out.states.push_back(source);
        return out;
    }

    std::queue<std::string> q;
    std::unordered_map<std::string, PathParent> parent;
    std::unordered_map<std::string, VectorRangeTreeMap> states_by_key;

    parent[src_key] = {"", "", 0};
    states_by_key[src_key] = source;
    q.push(src_key);

    while (!q.empty()) {
        std::string cur_key = q.front();
        q.pop();

        const PathParent cur_parent = parent[cur_key];
        const int cur_depth = cur_parent.depth;
        if (cur_depth >= depth_cap) {
            continue;
        }

        const VectorRangeTreeMap &cur = states_by_key[cur_key];
        std::vector<int> nodes(cur.original_nodes.begin(), cur.original_nodes.end());
        std::sort(nodes.begin(), nodes.end());

        for (int v : nodes) {
            if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                VectorRangeTreeMap nxt = cur;
                nxt.rotateLeft(v);
                std::string nxt_key = treeToString(nxt);
                if (!parent.count(nxt_key)) {
                    parent[nxt_key] = {cur_key, "rotateLeft(" + std::to_string(v) + ")", cur_depth + 1};
                    states_by_key[nxt_key] = std::move(nxt);
                    if (nxt_key == dst_key) {
                        out.found = true;
                        q = {};
                        break;
                    }
                    if (parent.size() >= state_cap) {
                        out.capped = true;
                        q = {};
                        break;
                    }
                    q.push(nxt_key);
                }
            }

            if (cur.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                VectorRangeTreeMap nxt = cur;
                nxt.rotateRight(v);
                std::string nxt_key = treeToString(nxt);
                if (!parent.count(nxt_key)) {
                    parent[nxt_key] = {cur_key, "rotateRight(" + std::to_string(v) + ")", cur_depth + 1};
                    states_by_key[nxt_key] = std::move(nxt);
                    if (nxt_key == dst_key) {
                        out.found = true;
                        q = {};
                        break;
                    }
                    if (parent.size() >= state_cap) {
                        out.capped = true;
                        q = {};
                        break;
                    }
                    q.push(nxt_key);
                }
            }
        }
    }

    if (!out.found) {
        return out;
    }

    std::vector<std::string> key_path;
    for (std::string cur = dst_key; !cur.empty(); cur = parent[cur].parent_key) {
        key_path.push_back(cur);
    }
    std::reverse(key_path.begin(), key_path.end());

    out.states.reserve(key_path.size());
    out.moves.reserve(key_path.size() > 0 ? key_path.size() - 1 : 0);
    for (size_t i = 0; i < key_path.size(); i++) {
        out.states.push_back(states_by_key[key_path[i]]);
        if (i > 0) {
            out.moves.push_back(parent[key_path[i]].move);
        }
    }

    return out;
}

// Definition: Print path steps and optional ASCII trees for each step
// Parameters: direction: a->b or b->a; T1/T2: source and target trees; dist: expected shortest distance; ascii: print tree shape if true
// Returns: nothing
// Errors: prints a skip/cap message if path search is too large
static void emitRotationPath(const char *direction,
                             const VectorRangeTreeMap &T1,
                             const VectorRangeTreeMap &T2,
                             int dist,
                             bool ascii) {
    if (dist < 0) {
        std::cout << "[PATH] direction=" << direction << " skipped (distance not found)\n";
        return;
    }

    // Keep debug path generation bounded so normal large runs do not explode
    const int n = static_cast<int>(T1.original_nodes.size());
    if (n > 14 || dist > 14) {
        std::cout << "[PATH] direction=" << direction
                  << " skipped (n/dist too large for debug path: n=" << n
                  << ", dist=" << dist << ")\n";
        return;
    }

    PathBuildResult path = buildShortestRotationPathBfs(T1, T2, dist, 250000);
    if (!path.found) {
        std::cout << "[PATH] direction=" << direction
                  << " not found within depth=" << dist;
        if (path.capped) std::cout << " (state cap reached)";
        std::cout << "\n";
        return;
    }

    std::cout << "[PATH] direction=" << direction
              << " steps=" << path.moves.size() << "\n";
    for (size_t i = 0; i < path.states.size(); i++) {
        std::cout << "  [step " << i << "]";
        if (i > 0) {
            std::cout << " move=" << path.moves[i - 1];
        } else {
            std::cout << " move=start";
        }
        std::cout << "\n";
        if (ascii) {
            printAsciiTreeDiagram(path.states[i], "    ");
        } else {
            std::cout << "    " << currentTreeToPreInString(path.states[i]) << "\n";
        }
        std::cout << "    relations:\n";
        printNodeRelations(path.states[i], "      ");
    }

    const VectorRangeTreeMap &result_tree = path.states.back();
    const bool matches_target = TreesEqual(result_tree, T2);
    std::cout << "[RESULT_TREE] direction=" << direction
              << " matches_target=" << (matches_target ? "true" : "false") << "\n";
    std::cout << "  resulting_pre_in=" << currentTreeToPreInString(result_tree) << "\n";
    std::cout << "  target_pre_in=" << currentTreeToPreInString(T2) << "\n";
    std::cout << "  resulting_struct=" << treeToString(result_tree) << "\n";
    std::cout << "  target_struct=" << treeToString(T2) << "\n";
    if (ascii) {
        std::cout << "  resulting_ascii:\n";
        printAsciiTreeDiagram(result_tree, "    ");
        std::cout << "  resulting_relations:\n";
        printNodeRelations(result_tree, "    ");
        std::cout << "  target_ascii:\n";
        printAsciiTreeDiagram(T2, "    ");
        std::cout << "  target_relations:\n";
        printNodeRelations(T2, "    ");
    }
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
            if (opts.print_trees) {
                std::cout << "[TREES] direction=" << direction << "\n";
                std::cout << "  source_pre_in=" << treeToPreInString(T1) << "\n";
                std::cout << "  target_pre_in=" << treeToPreInString(T2) << "\n";
                std::cout << "  source_struct=" << treeToString(T1) << "\n";
                std::cout << "  target_struct=" << treeToString(T2) << "\n";
            }

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

            if (opts.emit_path || opts.path_ascii) {
                const bool show_ascii = opts.path_ascii || opts.emit_path;
                emitRotationPath(direction, T1, T2, dist, show_ascii);
            }

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
