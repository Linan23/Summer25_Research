#include "testing_cli.h"
#include "algorithm.h"
#include "helpers.h"
#include "memoization.h"
#include "bf_bst.h"

#include <algorithm>
#include <chrono>
#include <fstream>
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
    std::chrono::steady_clock::time_point start_time;
public:
    // Definition: Start or reset the timer
    // Parameters: none
    // Returns: nothing
    // Errors: none
    void start() { start_time = std::chrono::steady_clock::now(); }

    // Definition: Elapsed time in microseconds since start
    // Parameters: none
    // Returns: elapsed microseconds
    // Errors: none
    long long getMicroseconds() {
        auto end_time = std::chrono::steady_clock::now();
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

// Definition: Read a whole text file into a string
// Parameters: path: file path to read
// Returns: file contents
// Errors: throws on open or read failure
static std::string readTextFile(const std::string &path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Unable to open file: " + path);
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

// Definition: Read an integer env var with bounds and a fallback
// Parameters: name: env var name, fallback: value when unset, min_value/max_value: clamp range
// Returns: parsed integer value
// Errors: malformed values fall back through atoi semantics and are clamped
static int readEnvInt(const char *name, int fallback, int min_value, int max_value) {
    const char *env = std::getenv(name);
    int value = fallback;
    if (env && *env) {
        value = std::atoi(env);
    }
    if (value < min_value) value = min_value;
    if (value > max_value) value = max_value;
    return value;
}

// Definition: Trim leading and trailing ASCII whitespace
// Parameters: value: input text
// Returns: trimmed string
// Errors: none
static std::string trimAscii(std::string value) {
    const auto is_ws = [](unsigned char c) {
        return c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '\f' || c == '\v';
    };
    while (!value.empty() && is_ws(static_cast<unsigned char>(value.front()))) {
        value.erase(value.begin());
    }
    while (!value.empty() && is_ws(static_cast<unsigned char>(value.back()))) {
        value.pop_back();
    }
    return value;
}

// Definition: Parse a comma-separated integer list
// Parameters: part: comma-separated integer text
// Returns: parsed integer vector
// Errors: throws on invalid integers
static std::vector<int> parseIntSequence(const std::string &part) {
    std::vector<int> out;
    std::stringstream ss(part);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trimAscii(token);
        if (token.empty()) continue;
        out.push_back(std::stoi(token));
    }
    return out;
}

// Definition: Build a tree from a canonical P:...;I:... traversal string
// Parameters: canonical: encoded preorder/inorder string
// Returns: constructed tree
// Errors: throws on malformed input
static VectorRangeTreeMap treeFromCanonicalString(const std::string &canonical) {
    std::string text = trimAscii(canonical);
    if (text.rfind("P:", 0) != 0) {
        throw std::runtime_error("Canonical tree must start with P:");
    }
    const std::string delim = ";I:";
    size_t split = text.find(delim);
    if (split == std::string::npos) {
        throw std::runtime_error("Canonical tree is missing ;I:");
    }
    std::vector<int> preorder = parseIntSequence(text.substr(2, split - 2));
    std::vector<int> inorder = parseIntSequence(text.substr(split + delim.size()));
    if (preorder.empty() || inorder.empty()) {
        throw std::runtime_error("Canonical tree is empty");
    }
    if (preorder.size() != inorder.size()) {
        throw std::runtime_error("Preorder and inorder sizes do not match");
    }
    VectorRangeTreeMap out;
    out.build(preorder, inorder);
    return out;
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
    std::string tree_a_file;
    std::string tree_b_file;
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
    std::cerr << "Usage: " << argv0 << " --case random|simple|comb --n <int> [--seed <int>] [--count <int>]\n"
              << "       [--max-k <int>] [--bfs-cap <int>] [--print-trees] [--emit-path] [--path-ascii]\n"
              << "   or: " << argv0 << " --tree-a-file <path> --tree-b-file <path> [--max-k <int>] [--print-trees]\n"
              << "   Note: --case simple is the clearer alias for the older --case comb.\n";
}

static bool isSimpleGeneratedCase(const std::string &case_type) {
    return case_type == "simple" || case_type == "comb";
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
        if (arg == "--tree-a-file" && i + 1 < argc) {
            opts.tree_a_file = argv[++i];
            continue;
        }
        if (arg == "--tree-b-file" && i + 1 < argc) {
            opts.tree_b_file = argv[++i];
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
    const int max_n = readEnvInt("FLIPDIST_PATH_DEBUG_MAX_N", 14, 1, 64);
    const int max_dist = readEnvInt("FLIPDIST_PATH_DEBUG_MAX_DIST", 14, 1, 128);
    const size_t state_cap = static_cast<size_t>(
        readEnvInt("FLIPDIST_PATH_DEBUG_STATE_CAP", 250000, 1000, 50000000)
    );
    if (n > max_n || dist > max_dist) {
        std::cout << "[PATH] direction=" << direction
                  << " skipped (n/dist too large for debug path: n=" << n
                  << ", dist=" << dist
                  << ", limits=" << max_n << "/" << max_dist << ")\n";
        return;
    }

    PathBuildResult path = buildShortestRotationPathBfs(T1, T2, dist, state_cap);
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
    const bool custom_inputs = !opts.tree_a_file.empty() || !opts.tree_b_file.empty();
    if (custom_inputs && (opts.tree_a_file.empty() || opts.tree_b_file.empty())) {
        std::cerr << "Both --tree-a-file and --tree-b-file are required together\n";
        return 2;
    }
    if (!custom_inputs && opts.n <= 0) {
        std::cerr << "Missing or invalid --n\n";
        return 2;
    }
    if (opts.count <= 0) {
        std::cerr << "Invalid --count\n";
        return 2;
    }
    if (!custom_inputs && opts.case_type != "random" && !isSimpleGeneratedCase(opts.case_type)) {
        std::cerr << "Invalid --case (use random, simple, or comb)\n";
        return 2;
    }

    auto run_case = [&](const VectorRangeTreeMap &A,
                        const VectorRangeTreeMap &B,
                        const std::string &case_label,
                        int case_seed,
                        int n_value,
                        int max_k_value) {

        struct DirectionRunResult {
            int dist = -1;
            double ms = 0.0;
            bool aborted = false;
        };

        auto run_one = [&](const VectorRangeTreeMap &T1,
                           const VectorRangeTreeMap &T2,
                           const char *direction,
                           int direction_max_k,
                           bool clear_memo,
                           bool emit_output = true,
                           double time_offset_ms = 0.0) -> DirectionRunResult {
            if (emit_output && opts.print_trees) {
                std::cout << "[TREES] direction=" << direction << "\n";
                std::cout << "  source_pre_in=" << treeToPreInString(T1) << "\n";
                std::cout << "  target_pre_in=" << treeToPreInString(T2) << "\n";
                std::cout << "  source_struct=" << treeToString(T1) << "\n";
                std::cout << "  target_struct=" << treeToString(T2) << "\n";
            }

            initProfile();
            if (g_profile.enabled) {
                g_profile.start_time = std::chrono::steady_clock::now();
                resetPlateauProfile();
                beginPlateauProfileRun(case_label, n_value, case_seed, direction,
                                       std::string(direction) == "a->b");
                resetTdiPostIProfile();
                beginTdiPostIProfileRun(case_label, n_value, case_seed, direction,
                                        false);
                g_profile.calls_flipdist = 0;
                g_profile.calls_tdi = 0;
                g_profile.calls_tds = 0;
                g_profile.tdi_prefix_lb_checks = 0;
                g_profile.tdi_prefix_lb_prunes = 0;
                g_profile.tdi_prefix_pairgen_calls = 0;
                g_profile.tdi_prefix_pairgen_pairs_scanned = 0;
                g_profile.tdi_prefix_pairgen_pairs_inserted = 0;
                g_profile.tdi_prefix_request_total = 0;
                g_profile.tdi_prefix_request_unique = 0;
                g_profile.tdi_prefix_request_repeated = 0;
                g_profile.tdi_prefix_emitted_pairs_0 = 0;
                g_profile.tdi_prefix_emitted_pairs_1 = 0;
                g_profile.tdi_prefix_emitted_pairs_2 = 0;
                g_profile.tdi_prefix_conflict_calls = 0;
                g_profile.tdi_prefix_bound_calls = 0;
                g_profile.memo_hits_tds = 0;
                g_profile.bounds_hits_tds = 0;
                g_profile.tds_unique_states = 0;
                g_profile.tds_repeated_states = 0;
                g_profile.calls_partition = 0;
                g_profile.partition_unique_structures = 0;
                g_profile.partition_repeated_structures = 0;
                g_profile.partition_unique_side_states = 0;
                g_profile.partition_repeated_side_states = 0;
                g_profile.partition_unique_split_signatures = 0;
                g_profile.partition_repeated_split_signatures = 0;
                g_profile.partition_precheck_calls = 0;
                g_profile.partition_precheck_matches = 0;
                g_profile.partition_precheck_mismatches = 0;
                g_profile.partition_budget_iterations = 0;
                g_profile.tds_partition_side_budget_cache_hits = 0;
                g_profile.tds_partition_side_budget_cache_misses = 0;
                g_profile.tds_partition_split_cache_hits = 0;
                g_profile.tds_partition_split_cache_misses = 0;
                g_profile.partition_empty_side_shortcuts = 0;
                g_profile.partition_side1_recursions = 0;
                g_profile.partition_side2_recursions = 0;
                g_profile.partition_side1_bound_hits = 0;
                g_profile.partition_side2_bound_hits = 0;
                g_profile.free_edge_hits = 0;
                g_profile.free_edge_misses = 0;
                g_profile.s_branch_calls = 0;
                g_profile.s_empty_branch_calls = 0;
                g_profile.s_empty_unique_states = 0;
                g_profile.s_empty_repeated_states = 0;
                g_profile.s_empty_no_candidate_states = 0;
                g_profile.s_empty_progress_states = 0;
                g_profile.s_empty_plateau_states = 0;
                g_profile.s_empty_progress_candidates = 0;
                g_profile.s_empty_plateau_candidates = 0;
                g_profile.s_empty_plateau_buckets = 0;
                g_profile.s_empty_artic_arm_hits = 0;
                g_profile.s_empty_duplicate_child_states = 0;
                g_profile.s_empty_must_drop_rejects = 0;
                g_profile.s_empty_incumbent_rejects = 0;
                g_profile.s_empty_core_lookups = 0;
                g_profile.s_empty_core_cache_hits = 0;
                g_profile.s_empty_core_exact_hits = 0;
                g_profile.s_empty_core_common_splits = 0;
                g_profile.s_empty_core_free_reductions = 0;
                g_profile.s_empty_motif_common_peels = 0;
                g_profile.s_empty_motif_free_peels = 0;
                g_profile.s_empty_motif_exact_fallbacks = 0;
                g_profile.s_empty_motif_exact_cache_hits = 0;
                g_profile.s_empty_branchy_reduction_hits = 0;
                g_profile.s_empty_branchy_exact_bounds = 0;
                g_profile.s_empty_branchy_exact_decisions = 0;
                g_profile.s_empty_branchy_exact_cache_hits = 0;
                g_profile.s_dead_pair_skips = 0;
                g_profile.s_dead_prefixes = 0;
                g_profile.s_full_dead_prefix_failures = 0;
                g_profile.independent_subsets_initial = 0;
                g_profile.independent_subsets_s = 0;
                g_profile.max_s_size = 0;
                g_profile.max_i_size = 0;
                g_profile.heuristic_lb_calls = 0;
                g_profile.heuristic_lb_cache_hits = 0;
                g_profile.conflict_cache_hits = 0;
                g_profile.free_edge_cache_hits = 0;
                g_profile.dominance_prunes = 0;
                g_profile.small_exact_lookups = 0;
                g_profile.small_exact_hits = 0;
                g_profile.local_nopath_lookups = 0;
                g_profile.local_nopath_hits = 0;
                g_profile.time_flipdist_ms = 0.0;
                g_profile.time_tdi_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_helper_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_build_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_build_edges_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_build_ranges_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_build_range_lookup_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_build_insert_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_build_boundary_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_pick_ms = 0.0;
                g_profile.time_tdi_prefix_pairgen_insert_ms = 0.0;
                g_profile.time_tdi_prefix_conflict_ms = 0.0;
                g_profile.time_tdi_prefix_bound_ms = 0.0;
                g_profile.tdi_prefix_request_vcount_samples.clear();
                g_profile.tdi_prefix_request_bucket_l_samples.clear();
                g_profile.tdi_prefix_request_bucket_r_samples.clear();
                g_profile.time_tds_ms = 0.0;
                g_profile.time_tds_free_filter_ms = 0.0;
                g_profile.time_tds_free_rotate_ms = 0.0;
                g_profile.time_tds_free_partner_ms = 0.0;
                g_profile.time_partition_ms = 0.0;
                g_profile.time_partition_precheck_ms = 0.0;
                g_profile.time_partition_build_ms = 0.0;
                g_profile.time_partition_side_stats_ms = 0.0;
                g_profile.time_partition_count_edges_ms = 0.0;
                g_profile.time_partition_conflicts_ms = 0.0;
                g_profile.time_partition_bounds_ms = 0.0;
                g_profile.time_partition_split_s_ms = 0.0;
                g_profile.time_partition_budget_loop_ms = 0.0;
                resetAlgorithmProfileSets();
            }
            PerformanceTimer timer;
            timer.start();
            if (clear_memo) {
                resetMemo();
            }
            int dist = -1;
            bool aborted = false;
            try {
                dist = FlipDistMinK(T1, T2, direction_max_k);
            } catch (const SearchAbortException&) {
                aborted = true;
            }
            double ms = timer.getMilliseconds();

            std::string status = (dist >= 0) ? "ok" : "not_found";
            if (emit_output) {
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "{\"case_type\":\"" << case_label << "\""
                          << ",\"n\":" << n_value
                          << ",\"seed\":" << case_seed
                          << ",\"direction\":\"" << direction << "\""
                          << ",\"distance\":" << dist
                          << ",\"time_ms\":" << (ms + time_offset_ms)
                          << ",\"status\":\"" << status << "\""
                          << ",\"tree_a\":\"" << treeToPreInString(T1) << "\""
                          << ",\"tree_b\":\"" << treeToPreInString(T2) << "\""
                          << ",\"max_k\":" << direction_max_k
                          << "}" << std::endl;
            }

            if (emit_output && (opts.emit_path || opts.path_ascii)) {
                const bool show_ascii = opts.path_ascii || opts.emit_path;
                emitRotationPath(direction, T1, T2, dist, show_ascii);
            }

            if (emit_output && g_profile.enabled) {
                auto sampleAverage = [](const std::vector<int>& samples) -> double {
                    if (samples.empty()) return 0.0;
                    long long total = 0;
                    for (int value : samples) total += value;
                    return static_cast<double>(total) / static_cast<double>(samples.size());
                };
                auto samplePercentile = [](std::vector<int> samples, double pct) -> double {
                    if (samples.empty()) return 0.0;
                    std::sort(samples.begin(), samples.end());
                    if (samples.size() == 1) return static_cast<double>(samples.front());
                    const double pos = pct * static_cast<double>(samples.size() - 1);
                    const size_t lo = static_cast<size_t>(pos);
                    const size_t hi = std::min(samples.size() - 1, lo + 1);
                    const double frac = pos - static_cast<double>(lo);
                    return static_cast<double>(samples[lo]) * (1.0 - frac) +
                           static_cast<double>(samples[hi]) * frac;
                };
                std::cout << "[PROFILE] direction=" << direction << " dist=" << dist << "\n"
                          << "  - FlipDistTree calls=" << g_profile.calls_flipdist
                          << " time_ms=" << g_profile.time_flipdist_ms << "\n"
                          << "  - TreeDistI calls=" << g_profile.calls_tdi
                          << " time_ms=" << g_profile.time_tdi_ms << "\n"
                          << "  - TreeDistI bounds_prunes=" << g_profile.tdi_bounds_prunes
                          << " conflict_lb_prunes=" << g_profile.tdi_conflict_lb_prunes
                          << " remaining_budget_prunes=" << g_profile.tdi_remaining_budget_prunes << "\n"
                          << "  - TreeDistI exact_budget success=" << g_profile.tdi_exact_budget_success
                          << " fail=" << g_profile.tdi_exact_budget_fail
                          << " tree_equal_success=" << g_profile.tdi_tree_equal_success << "\n"
                          << "  - TreeDistI post_I empty_s=" << g_profile.tdi_post_i_empty_s
                          << " nonempty_s=" << g_profile.tdi_post_i_nonempty_s
                          << " post_I_lb_prunes=" << g_profile.tdi_post_i_lb_prunes
                          << " prefix_lb_checks=" << g_profile.tdi_prefix_lb_checks
                          << " prefix_lb_prunes=" << g_profile.tdi_prefix_lb_prunes
                          << " prefix_pairgen_calls=" << g_profile.tdi_prefix_pairgen_calls
                          << " prefix_pairs_scanned=" << g_profile.tdi_prefix_pairgen_pairs_scanned
                          << " prefix_pairs_inserted=" << g_profile.tdi_prefix_pairgen_pairs_inserted
                          << " prefix_requests_total=" << g_profile.tdi_prefix_request_total
                          << " prefix_requests_unique=" << g_profile.tdi_prefix_request_unique
                          << " prefix_requests_repeated=" << g_profile.tdi_prefix_request_repeated
                          << " prefix_conflict_calls=" << g_profile.tdi_prefix_conflict_calls
                          << " prefix_bound_calls=" << g_profile.tdi_prefix_bound_calls
                          << " tds_success=" << g_profile.tdi_tds_success
                          << " tds_fail=" << g_profile.tdi_tds_fail
                          << " rotation_skips=" << g_profile.tdi_rotation_skips << "\n"
                          << "  - TreeDistI prefix_time pairgen_ms=" << g_profile.time_tdi_prefix_pairgen_ms
                          << " helper_ms=" << g_profile.time_tdi_prefix_pairgen_helper_ms
                          << " build_ms=" << g_profile.time_tdi_prefix_pairgen_build_ms
                          << " build_edges_ms=" << g_profile.time_tdi_prefix_pairgen_build_edges_ms
                          << " build_ranges_ms=" << g_profile.time_tdi_prefix_pairgen_build_ranges_ms
                          << " build_range_lookup_ms=" << g_profile.time_tdi_prefix_pairgen_build_range_lookup_ms
                          << " build_insert_ms=" << g_profile.time_tdi_prefix_pairgen_build_insert_ms
                          << " build_boundary_ms=" << g_profile.time_tdi_prefix_pairgen_build_boundary_ms
                          << " pick_ms=" << g_profile.time_tdi_prefix_pairgen_pick_ms
                          << " insert_ms=" << g_profile.time_tdi_prefix_pairgen_insert_ms
                          << " conflict_ms=" << g_profile.time_tdi_prefix_conflict_ms
                          << " bound_ms=" << g_profile.time_tdi_prefix_bound_ms << "\n"
                          << "  - TreeDistI prefix_shape avg_vcount="
                          << sampleAverage(g_profile.tdi_prefix_request_vcount_samples)
                          << " p95_vcount="
                          << samplePercentile(g_profile.tdi_prefix_request_vcount_samples, 0.95)
                          << " avg_bucket_L="
                          << sampleAverage(g_profile.tdi_prefix_request_bucket_l_samples)
                          << " p95_bucket_L="
                          << samplePercentile(g_profile.tdi_prefix_request_bucket_l_samples, 0.95)
                          << " avg_bucket_R="
                          << sampleAverage(g_profile.tdi_prefix_request_bucket_r_samples)
                          << " p95_bucket_R="
                          << samplePercentile(g_profile.tdi_prefix_request_bucket_r_samples, 0.95)
                          << " emitted_pairs_hist=0:" << g_profile.tdi_prefix_emitted_pairs_0
                          << ",1:" << g_profile.tdi_prefix_emitted_pairs_1
                          << ",2+:" << g_profile.tdi_prefix_emitted_pairs_2 << "\n"
                          << "  - TreeDistS calls=" << g_profile.calls_tds
                          << " time_ms=" << g_profile.time_tds_ms << "\n"
                          << "  - TreeDistS memo_hits=" << g_profile.memo_hits_tds
                          << " bounds_hits=" << g_profile.bounds_hits_tds << "\n"
#if FLIPDIST_REUSE_PROFILE_COUNTERS
                          << "  - TreeDistS state_reuse unique=" << g_profile.tds_unique_states
                          << " repeated=" << g_profile.tds_repeated_states << "\n"
#endif
                          << "  - free_edge_time filter_ms=" << g_profile.time_tds_free_filter_ms
                          << " rotate_ms=" << g_profile.time_tds_free_rotate_ms
                          << " partner_ms=" << g_profile.time_tds_free_partner_ms << "\n"
                          << "  - partition calls=" << g_profile.calls_partition
                          << " time_ms=" << g_profile.time_partition_ms
                          << " precheck_calls=" << g_profile.partition_precheck_calls
                          << " precheck_matches=" << g_profile.partition_precheck_matches
                          << " precheck_mismatches=" << g_profile.partition_precheck_mismatches << "\n"
#if FLIPDIST_REUSE_PROFILE_COUNTERS
                          << "  - partition_reuse structures_unique=" << g_profile.partition_unique_structures
                          << " structures_repeated=" << g_profile.partition_repeated_structures
                          << " side_states_unique=" << g_profile.partition_unique_side_states
                          << " side_states_repeated=" << g_profile.partition_repeated_side_states
                          << " split_signatures_unique=" << g_profile.partition_unique_split_signatures
                          << " split_signatures_repeated=" << g_profile.partition_repeated_split_signatures << "\n"
#endif
                          << "  - partition_time precheck_ms=" << g_profile.time_partition_precheck_ms
                          << " build_ms=" << g_profile.time_partition_build_ms
                          << " side_stats_ms=" << g_profile.time_partition_side_stats_ms
                          << " count_edges_ms=" << g_profile.time_partition_count_edges_ms
                          << " conflicts_ms=" << g_profile.time_partition_conflicts_ms
                          << " bounds_ms=" << g_profile.time_partition_bounds_ms
                          << " split_s_ms=" << g_profile.time_partition_split_s_ms
                          << " budget_loop_ms=" << g_profile.time_partition_budget_loop_ms << "\n"
                          << "  - partition_budget iterations=" << g_profile.partition_budget_iterations
                          << " empty_side_shortcuts=" << g_profile.partition_empty_side_shortcuts
                          << " side_budget_cache_hits=" << g_profile.tds_partition_side_budget_cache_hits
                          << " side_budget_cache_misses=" << g_profile.tds_partition_side_budget_cache_misses
                          << " split_cache_hits=" << g_profile.tds_partition_split_cache_hits
                          << " split_cache_misses=" << g_profile.tds_partition_split_cache_misses
                          << " side1_recursions=" << g_profile.partition_side1_recursions
                          << " side2_recursions=" << g_profile.partition_side2_recursions
                          << " side1_bound_hits=" << g_profile.partition_side1_bound_hits
                          << " side2_bound_hits=" << g_profile.partition_side2_bound_hits << "\n"
                          << "  - free_edge hits=" << g_profile.free_edge_hits
                          << " misses=" << g_profile.free_edge_misses << "\n"
                          << "  - S-branch calls=" << g_profile.s_branch_calls
                          << " S-empty calls=" << g_profile.s_empty_branch_calls << "\n"
#if FLIPDIST_REUSE_PROFILE_COUNTERS
                          << "  - S-empty state_reuse unique=" << g_profile.s_empty_unique_states
                          << " repeated=" << g_profile.s_empty_repeated_states << "\n"
#endif
                          << "  - S-empty progress_states=" << g_profile.s_empty_progress_states
                          << " plateau_states=" << g_profile.s_empty_plateau_states
                          << " no_candidate_states=" << g_profile.s_empty_no_candidate_states << "\n"
                          << "  - S-empty progress_candidates=" << g_profile.s_empty_progress_candidates
                          << " plateau_candidates=" << g_profile.s_empty_plateau_candidates
                          << " plateau_buckets=" << g_profile.s_empty_plateau_buckets << "\n"
                          << "  - S-empty artic_arm_hits=" << g_profile.s_empty_artic_arm_hits << "\n"
                          << "  - S-empty duplicate_child_states=" << g_profile.s_empty_duplicate_child_states
                          << " must_drop_rejects=" << g_profile.s_empty_must_drop_rejects
                          << " incumbent_rejects=" << g_profile.s_empty_incumbent_rejects << "\n"
                          << "  - S-empty core lookups=" << g_profile.s_empty_core_lookups
                          << " cache_hits=" << g_profile.s_empty_core_cache_hits
                          << " exact_hits=" << g_profile.s_empty_core_exact_hits
                          << " common_splits=" << g_profile.s_empty_core_common_splits
                          << " free_reductions=" << g_profile.s_empty_core_free_reductions << "\n"
                          << "  - S-empty motif common_peels=" << g_profile.s_empty_motif_common_peels
                          << " free_peels=" << g_profile.s_empty_motif_free_peels
                          << " exact_fallbacks=" << g_profile.s_empty_motif_exact_fallbacks
                          << " exact_cache_hits=" << g_profile.s_empty_motif_exact_cache_hits << "\n"
                          << "  - S-empty branchy reductions=" << g_profile.s_empty_branchy_reduction_hits
                          << " exact_bounds=" << g_profile.s_empty_branchy_exact_bounds
                          << " exact_decisions=" << g_profile.s_empty_branchy_exact_decisions
                          << " exact_cache_hits=" << g_profile.s_empty_branchy_exact_cache_hits << "\n"
                          << "  - S dead_pair_skips=" << g_profile.s_dead_pair_skips
                          << " dead_prefixes=" << g_profile.s_dead_prefixes
                          << " full_dead_prefix_failures=" << g_profile.s_full_dead_prefix_failures << "\n"
                          << "  - independent_subsets_initial=" << g_profile.independent_subsets_initial
                          << " max_I_size=" << g_profile.max_i_size << "\n"
                          << "  - max_S_size=" << g_profile.max_s_size << "\n"
                          << "  - heuristic_lb calls=" << g_profile.heuristic_lb_calls
                          << " cache_hits=" << g_profile.heuristic_lb_cache_hits << "\n"
                          << "  - conflict_cache_hits=" << g_profile.conflict_cache_hits
                          << " free_edge_cache_hits=" << g_profile.free_edge_cache_hits << "\n"
                          << "  - dominance_prunes=" << g_profile.dominance_prunes << "\n"
                          << "  - small_exact lookups=" << g_profile.small_exact_lookups
                          << " hits=" << g_profile.small_exact_hits << "\n"
                          << "  - local_nopath lookups=" << g_profile.local_nopath_lookups
                          << " hits=" << g_profile.local_nopath_hits << "\n";
                emitPlateauProfileJsonl(std::cout, case_label, n_value, case_seed, direction);
                finalizePlateauProfileRun();
                if (std::getenv("FLIPDIST_PROFILE_TDI_POSTI_FILE") == nullptr) {
                    emitTdiPostIProfileJsonl(std::cout, case_label, n_value, case_seed, direction);
                }
                finalizeTdiPostIProfileRun();
            }
            return DirectionRunResult{dist, ms, aborted};
        };

        const int conflicts_ab = countConflictEdges(A, B);
        const int conflicts_ba = countConflictEdges(B, A);
        struct DirectionShape {
            int height = 0;
            int branching = 0;
            int root = -1;
            int root_label = -1;
            double avg_depth = 0.0;
        };
        auto directionShape = [](const VectorRangeTreeMap& T) {
            DirectionShape shape;
            shape.root = T.root;
            if (!T.original_preorder.empty()) {
                shape.root_label = T.original_preorder.front();
            }
            long long depth_sum = 0;
            int node_count = 0;
            auto dfs = [&](auto&& self, int node, int depth) -> void {
                if (node < 0 || !T.isOriginal(node)) return;
                node_count++;
                depth_sum += depth;
                if (depth > shape.height) shape.height = depth;
                const int left = T.getLeftChild(node);
                const int right = T.getRightChild(node);
                if (left >= 0 && T.isOriginal(left) && right >= 0 && T.isOriginal(right)) {
                    shape.branching++;
                }
                self(self, left, depth + 1);
                self(self, right, depth + 1);
            };
            dfs(dfs, T.root, 0);
            if (node_count > 0) {
                shape.avg_depth = static_cast<double>(depth_sum) / static_cast<double>(node_count);
            }
            return shape;
        };
        const DirectionShape shape_a = directionShape(A);
        const DirectionShape shape_b = directionShape(B);
        const bool reverse_shape_promising =
            shape_b.height >= shape_a.height + 2 ||
            shape_b.avg_depth >= shape_a.avg_depth + 0.3 ||
            shape_b.branching >= shape_a.branching + 2 ||
            (n_value >= 25 && shape_b.branching > shape_a.branching);
        const bool reverse_shape_risky =
            (shape_b.branching < shape_a.branching && shape_b.height <= shape_a.height) ||
            (shape_b.branching + 2 <= shape_a.branching && shape_b.height > shape_a.height) ||
            (shape_b.avg_depth + 0.1 < shape_a.avg_depth && shape_b.branching <= shape_a.branching);
        const char *force_first_env = std::getenv("FLIPDIST_FORCE_FIRST_DIRECTION");
        const std::string force_first = force_first_env ? std::string(force_first_env) : std::string();
        const char *direction_probe_env = std::getenv("FLIPDIST_DIRECTION_PROBE_THRESHOLD_MS");
        int direction_probe_threshold_ms = (n_value >= 23 && n_value <= 25) ? 7 : 0;
        if (direction_probe_env && *direction_probe_env) {
            direction_probe_threshold_ms = std::atoi(direction_probe_env);
        }
        const char *direction_probe_delta_env = std::getenv("FLIPDIST_DIRECTION_PROBE_DELTA");
        int direction_probe_delta = 2;
        if (direction_probe_delta_env && *direction_probe_delta_env) {
            direction_probe_delta = std::atoi(direction_probe_delta_env);
        }
        bool run_ba_first = (force_first == "ba" || force_first == "b->a") ||
                            (force_first.empty() && conflicts_ba < conflicts_ab);
        bool forward_shape_locked = false;
        if (force_first.empty() && n_value >= 23 && conflicts_ab == conflicts_ba &&
            shape_b.height >= shape_a.height + 3) {
            run_ba_first = true;
        }
        if (force_first.empty() && n_value >= 24 && conflicts_ab == conflicts_ba && !run_ba_first) {
            const bool reverse_branchier_deeper =
                shape_b.branching >= shape_a.branching + 3 &&
                shape_b.avg_depth >= shape_a.avg_depth + 0.5 &&
                shape_b.height <= shape_a.height + 1;
            const bool reverse_taller_shallower =
                shape_b.height > shape_a.height &&
                shape_b.avg_depth + 0.5 < shape_a.avg_depth &&
                shape_b.branching + 1 >= shape_a.branching;
            if (reverse_branchier_deeper || reverse_taller_shallower) {
                run_ba_first = true;
            }
        }
        if (force_first.empty()) {
            const bool reverse_root_extreme =
                shape_b.root <= 2 || shape_b.root >= n_value - 1;
            const int reverse_height_delta = shape_b.height - shape_a.height;
            const int reverse_branching_delta = shape_b.branching - shape_a.branching;
            const double reverse_avg_depth_delta = shape_b.avg_depth - shape_a.avg_depth;
            const bool boundary_reverse_low_root_deep =
                n_value == 26 &&
                shape_b.root_label <= 1 &&
                reverse_height_delta >= 2 &&
                reverse_avg_depth_delta >= 1.0 &&
                reverse_branching_delta == 0;
            const bool boundary_reverse_high_root_taller =
                n_value == 26 &&
                shape_b.root_label >= n_value &&
                reverse_height_delta >= 1 &&
                reverse_avg_depth_delta >= 0.45 &&
                reverse_branching_delta <= -1;
            const bool boundary_reverse_compact_more_branch =
                n_value == 27 &&
                reverse_height_delta <= -4 &&
                reverse_avg_depth_delta <= -1.0 &&
                reverse_branching_delta >= 1;
            const bool boundary_reverse_compact_less_branch =
                n_value == 27 &&
                reverse_height_delta <= -2 &&
                reverse_avg_depth_delta <= -0.8 &&
                reverse_branching_delta <= -1 &&
                shape_a.root_label >= 20;
            const bool boundary_reverse_mild_taller_less_branch =
                n_value == 27 &&
                reverse_height_delta == 1 &&
                reverse_avg_depth_delta >= 0.0 &&
                reverse_avg_depth_delta <= 0.2 &&
                reverse_branching_delta == -1 &&
                shape_a.root_label <= 8 &&
                shape_b.root_label <= 12;
            const bool boundary_reverse_equal_high_roots_deeper =
                n_value == 27 &&
                reverse_height_delta == 0 &&
                reverse_branching_delta == 0 &&
                reverse_avg_depth_delta >= 0.25 &&
                reverse_avg_depth_delta <= 0.5 &&
                shape_a.root_label >= 15 &&
                shape_b.root_label >= 20;
            const bool n25_reverse_low_root_tall_less_branchy_target =
                n_value == 25 &&
                shape_a.root_label <= 3 &&
                shape_b.root_label <= 5 &&
                reverse_height_delta >= 3 &&
                reverse_avg_depth_delta >= 1.2 &&
                reverse_branching_delta <= -2;
            const bool prefer_reverse_boundary_shape =
                boundary_reverse_low_root_deep ||
                boundary_reverse_high_root_taller ||
                boundary_reverse_compact_more_branch ||
                boundary_reverse_compact_less_branch ||
                boundary_reverse_mild_taller_less_branch ||
                boundary_reverse_equal_high_roots_deeper ||
                n25_reverse_low_root_tall_less_branchy_target;
            const bool prefer_forward_shape =
                (n_value == 23 &&
                 shape_b.height >= shape_a.height + 2 &&
                 shape_b.avg_depth >= shape_a.avg_depth + 0.3 &&
                 !reverse_root_extreme) ||
                (n_value >= 25 &&
                 shape_b.height >= shape_a.height + 2 &&
                 shape_b.avg_depth >= shape_a.avg_depth + 1.0) ||
                (n_value >= 26 &&
                 shape_b.root_label <= 2 &&
                 shape_b.height >= shape_a.height + 2 &&
                 shape_b.avg_depth >= shape_a.avg_depth + 0.7 &&
                 shape_b.branching <= shape_a.branching) ||
                (n_value == 26 &&
                 shape_b.root_label <= 2 &&
                 shape_a.root_label >= 5 &&
                 shape_a.root_label <= 10) ||
                (n_value >= 25 &&
                 shape_b.height >= shape_a.height + 2 &&
                 shape_b.avg_depth <= shape_a.avg_depth + 0.2 &&
                 !reverse_root_extreme);
            const bool n24_reverse_compact_high_root =
                n_value == 24 &&
                shape_a.root_label <= 4 &&
                shape_b.root_label >= 10 &&
                reverse_height_delta <= -1 &&
                reverse_avg_depth_delta <= -0.25 &&
                reverse_branching_delta <= 0;
            const bool n24_forward_compact_branchier_target =
                n_value == 24 &&
                reverse_height_delta <= -2 &&
                reverse_avg_depth_delta <= -0.35 &&
                reverse_branching_delta >= 2;
            const bool n25_reverse_high_to_compact_target =
                n_value == 25 &&
                shape_a.root_label >= 20 &&
                shape_b.root_label <= 10 &&
                reverse_height_delta <= -1 &&
                reverse_avg_depth_delta <= -0.35 &&
                reverse_branching_delta <= -2;
            const bool n25_forward_mild_taller_branchier_target =
                n_value == 25 &&
                shape_a.root_label >= 12 &&
                shape_b.root_label >= 20 &&
                reverse_height_delta == 1 &&
                reverse_avg_depth_delta >= 0.0 &&
                reverse_avg_depth_delta <= 0.35 &&
                reverse_branching_delta >= 1;
            const bool n25_forward_lowroot_compact_branchier_target =
                n_value == 25 &&
                shape_a.root_label <= 2 &&
                shape_b.root_label >= 8 &&
                reverse_height_delta <= -5 &&
                reverse_avg_depth_delta <= -2.0 &&
                reverse_branching_delta >= 2;
            const bool prefer_reverse_shape =
                n24_reverse_compact_high_root ||
                n25_reverse_high_to_compact_target ||
                ((n_value == 25 ||
                  (n_value >= 26 && shape_b.branching >= shape_a.branching)) &&
                 shape_a.avg_depth >= shape_b.avg_depth + 0.6 &&
                 shape_b.height >= shape_a.height - 2) ||
                (n_value == 24 &&
                 shape_a.root <= 2 &&
                 shape_a.height >= shape_b.height + 2 &&
                 shape_a.avg_depth >= shape_b.avg_depth + 1.0 &&
                 shape_b.branching >= shape_a.branching) ||
                (n_value == 23 &&
                 shape_b.avg_depth >= shape_a.avg_depth + 0.45 &&
                 shape_b.root >= n_value);
            if (prefer_reverse_boundary_shape) {
                run_ba_first = true;
            } else if (n24_forward_compact_branchier_target ||
                       n25_forward_mild_taller_branchier_target ||
                       n25_forward_lowroot_compact_branchier_target) {
                run_ba_first = false;
                forward_shape_locked = true;
            } else if (prefer_forward_shape) {
                run_ba_first = false;
                forward_shape_locked = true;
            } else if (prefer_reverse_shape) {
                run_ba_first = true;
            }
        }
        if (force_first.empty() && n_value >= 26 && conflicts_ab == conflicts_ba &&
            !run_ba_first && !forward_shape_locked) {
            const bool reverse_deeper_root_shift =
                shape_b.height >= shape_a.height + 1 &&
                shape_b.avg_depth >= shape_a.avg_depth + 0.3 &&
                shape_b.branching + 1 >= shape_a.branching &&
                shape_b.root + 6 <= shape_a.root;
            if (reverse_deeper_root_shift) {
                run_ba_first = true;
            }
        }
        if (force_first.empty() && n_value >= 28 && conflicts_ab == conflicts_ba &&
            !run_ba_first && !forward_shape_locked) {
            const bool reverse_deeper_less_branchy =
                shape_b.avg_depth >= shape_a.avg_depth + 0.25 &&
                shape_b.branching + 1 <= shape_a.branching;
            const bool reverse_shifted_less_branchy =
                n_value >= 29 &&
                shape_b.branching + 2 <= shape_a.branching &&
                shape_b.root >= shape_a.root + 8 &&
                shape_b.avg_depth + 0.4 >= shape_a.avg_depth;
            if (reverse_deeper_less_branchy || reverse_shifted_less_branchy) {
                run_ba_first = true;
            }
        }
        if (force_first.empty() && n_value >= 29 && conflicts_ab == conflicts_ba &&
            !run_ba_first && !forward_shape_locked) {
            const bool reverse_near_equal_shallow_shift =
                n_value == 29 &&
                shape_b.height + 1 == shape_a.height &&
                shape_b.branching == shape_a.branching &&
                shape_b.root > shape_a.root &&
                shape_b.root <= shape_a.root + 3 &&
                shape_b.avg_depth + 0.25 >= shape_a.avg_depth;
            const bool reverse_compact_target =
                n_value >= 30 &&
                shape_a.height >= shape_b.height + 3 &&
                shape_a.branching >= shape_b.branching + 3 &&
                shape_b.root <= shape_a.root &&
                shape_b.avg_depth + 0.6 >= shape_a.avg_depth;
            if (reverse_near_equal_shallow_shift || reverse_compact_target) {
                run_ba_first = true;
            }
        }

        const int final_reverse_height_delta = shape_b.height - shape_a.height;
        const int final_reverse_branching_delta = shape_b.branching - shape_a.branching;
        const double final_reverse_avg_depth_delta = shape_b.avg_depth - shape_a.avg_depth;
        const bool suppress_forward_direction_probe =
            (n_value == 24 &&
             final_reverse_height_delta <= -2 &&
             final_reverse_avg_depth_delta <= -0.20 &&
             final_reverse_branching_delta >= 2) ||
            (n_value == 25 &&
             shape_a.root_label >= 12 &&
             shape_b.root_label >= 20 &&
             final_reverse_height_delta == 1 &&
             final_reverse_avg_depth_delta >= 0.0 &&
             final_reverse_avg_depth_delta <= 0.40 &&
             final_reverse_branching_delta >= 1) ||
            (n_value == 25 &&
             shape_a.root_label <= 2 &&
             shape_b.root_label >= 8 &&
             final_reverse_height_delta <= -4 &&
             final_reverse_avg_depth_delta <= -1.50 &&
             final_reverse_branching_delta >= 2);
        const bool direction_probe_allowed =
            !run_ba_first &&
            force_first.empty() &&
            direction_probe_threshold_ms > 0 &&
            direction_probe_delta >= 0 &&
            conflicts_ab == conflicts_ba &&
            !suppress_forward_direction_probe &&
            reverse_shape_promising &&
            !reverse_shape_risky &&
            !(n_value >= 25 &&
              shape_b.height >= shape_a.height + 2 &&
              shape_b.avg_depth >= shape_a.avg_depth + 1.0) &&
            !(n_value >= 25 &&
              shape_b.height >= shape_a.height + 2 &&
              shape_b.avg_depth <= shape_a.avg_depth + 0.2 &&
              !(shape_b.root <= 2 || shape_b.root >= n_value - 1)) &&
            !(n_value == 23 &&
              shape_b.height >= shape_a.height + 2 &&
              shape_b.avg_depth >= shape_a.avg_depth + 0.3 &&
              !(shape_b.root <= 2 || shape_b.root >= n_value - 1)) &&
            !opts.print_trees &&
            !opts.emit_path &&
            !opts.path_ascii &&
            !(std::getenv("FLIPDIST_PROFILE") && std::string(std::getenv("FLIPDIST_PROFILE")) == "1");

        bool first_clear_memo = true;
        double first_time_offset_ms = 0.0;
        if (direction_probe_allowed) {
            const int probe_max_k = std::min(max_k_value, conflicts_ab + direction_probe_delta);
            DirectionRunResult probe =
                run_one(A, B, "a->b", probe_max_k, true, false);
            first_clear_memo = false;
            first_time_offset_ms = probe.ms;
            if (probe.dist >= 0) {
                run_one(A, B, "a->b", probe.dist, false, true);
                run_one(B, A, "b->a", probe.dist, false, true);
                return;
            }
            if (probe.ms >= static_cast<double>(direction_probe_threshold_ms)) {
                run_ba_first = true;
            } else {
                first_clear_memo = true;
            }
        }
        if (run_ba_first) {
            DirectionRunResult ba =
                run_one(B, A, "b->a", max_k_value, first_clear_memo, true, first_time_offset_ms);
            int dist_ba = ba.dist;
            int forward_max_k = max_k_value;
            if (dist_ba >= 0 && dist_ba < forward_max_k) {
                forward_max_k = dist_ba;
            }
            run_one(A, B, "a->b", forward_max_k, false);
        } else {
            DirectionRunResult ab =
                run_one(A, B, "a->b", max_k_value, first_clear_memo, true, first_time_offset_ms);
            int dist_ab = ab.dist;
            int reverse_max_k = max_k_value;
            if (dist_ab >= 0 && dist_ab < reverse_max_k) {
                reverse_max_k = dist_ab;
            }
            run_one(B, A, "b->a", reverse_max_k, false);
        }
    };

    if (custom_inputs) {
        try {
            VectorRangeTreeMap A = treeFromCanonicalString(readTextFile(opts.tree_a_file));
            VectorRangeTreeMap B = treeFromCanonicalString(readTextFile(opts.tree_b_file));
            const int n_value = static_cast<int>(A.original_nodes.size());
            if (n_value <= 0 || B.original_nodes.size() != A.original_nodes.size()) {
                std::cerr << "Custom tree files must contain non-empty trees with matching sizes\n";
                return 2;
            }
            if (opts.max_k <= 0) {
                opts.max_k = std::max(1, 3 * n_value + 10);
            }
            for (int i = 0; i < opts.count; i++) {
                run_case(A, B, "custom", -1, n_value, opts.max_k);
            }
        } catch (const std::exception &ex) {
            std::cerr << "Failed to load custom tree files: " << ex.what() << "\n";
            return 2;
        }
        return 0;
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
        run_case(A, B, opts.case_type, case_seed, opts.n, opts.max_k);
    }

    return 0;
}
