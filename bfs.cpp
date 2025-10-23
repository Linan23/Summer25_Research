// Centralises all breadth-first search variants used by the project. The file
// implements the plain FIFO solver, the heuristic-enhanced version, and the
// hashed bidirectional meet-in-the-middle search. The modern default flows via
// the optimised solver, automatically enabling rotation filters on large
// non-comb inputs and switching to BiBFS whenever time/space caps are hit.
// Helper utilities at the top of the file (canonical hashing, heuristics,
// env-flag readers) are shared across both C++ services and the comparison
// harness.

#include "rotation_tree.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>


// Returns a deterministic string encoding of a tree; helpful for diagnostics and
// collision checks because it orders nodes by their original label.
std::string canonicalKey(const VectorRangeTreeMap& T) {
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());

    std::string s;
    s.reserve(nodes.size() * 12);

    auto push32 = [&](int x){
        for (int k = 3; k >= 0; --k)
            s.push_back(static_cast<char>(static_cast<uint8_t>((x >> (k*8)) & 0xFF)));
    };

    s.push_back('R'); push32(T.root);
    s.push_back('|');

    for (int v : nodes) {
        push32(v);
        push32(T.getLeftChild(v));
        push32(T.getRightChild(v));
    }
    return s;
}

// Computes a 64-bit FNV-1a style hash that matches the information in
// canonicalKey but avoids heap allocations during the hot BFS loops.
uint64_t canonicalHash(const VectorRangeTreeMap& T) {
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());

    uint64_t h = 1469598103934665603ULL;
    auto mixByte = [&](uint8_t b){
        h ^= b;
        h *= 1099511628211ULL;
    };
    auto mixInt = [&](int value){
        uint32_t x = static_cast<uint32_t>(static_cast<int32_t>(value));
        mixByte(static_cast<uint8_t>(x & 0xFF));
        mixByte(static_cast<uint8_t>((x >> 8) & 0xFF));
        mixByte(static_cast<uint8_t>((x >> 16) & 0xFF));
        mixByte(static_cast<uint8_t>((x >> 24) & 0xFF));
    };

    mixByte(static_cast<uint8_t>('R'));
    mixInt(T.root);
    mixByte(static_cast<uint8_t>('|'));

    for (int v : nodes) {
        mixInt(v);
        mixInt(T.getLeftChild(v));
        mixInt(T.getRightChild(v));
    }

    return h;
}


namespace {
using EdgeSet = std::unordered_set<std::pair<int,int>,PairHash,PairEq>;

uint64_t mix64(uint64_t h, uint64_t k) {
    h ^= k + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    return h;
}

uint64_t subtreeHash(const VectorRangeTreeMap& tree, int node,
                     std::unordered_map<int,uint64_t>& memo)
{
    if (node == VectorRangeTreeMap::NO_CHILD || !tree.isOriginal(node))
        return 0x9e3779b97f4a7c15ULL;
    if (auto it = memo.find(node); it != memo.end()) return it->second;

    uint64_t h = 1469598103934665603ULL;
    h = mix64(h, static_cast<uint64_t>(node));

    int L = tree.getLeftChild(node);
    int R = tree.getRightChild(node);

    uint64_t hL = subtreeHash(tree, L, memo);
    uint64_t hR = subtreeHash(tree, R, memo);
    h = mix64(h, hL);
    h = mix64(h, hR);
    memo.emplace(node, h);
    return h;
}

// Materialises the parent->child edge set for a tree; used by the rotation filter
// to decide whether a move introduces a target edge.
EdgeSet collectEdgeSet(const VectorRangeTreeMap& T) {
    EdgeSet edges;
    T.collectEdges(T.root, edges);
    return edges;
}

// Returns true if rotating at node v will create at least one parent-child edge
// that belongs to the target tree according to target_edges.
bool rotationIntroducesTargetEdge(const VectorRangeTreeMap& cur,
                                  const EdgeSet& target_edges,
                                  int v,
                                  bool left_rotation)
{
    using TV = VectorRangeTreeMap;
    const int child = left_rotation ? cur.getRightChild(v)
                                    : cur.getLeftChild(v);
    if (child == TV::NO_CHILD || !cur.isOriginal(child)) return false;

    const int parent = cur.getParent(v);
    if (parent != TV::NO_PARENT && cur.isOriginal(parent)) {
        if (target_edges.count({parent, child})) return true;
    }

    if (target_edges.count({child, v})) return true;

    const int transferred = left_rotation ? cur.getLeftChild(child)
                                          : cur.getRightChild(child);
    if (transferred != TV::NO_CHILD && cur.isOriginal(transferred)) {
        if (target_edges.count({v, transferred})) return true;
    }

    return false;
}

// Counts how many parent->child edges cur already shares with the target.
int countMatchingTargetEdges(const VectorRangeTreeMap& tree,
                             const EdgeSet& target_edges)
{
    if (target_edges.empty()) return 0;
    EdgeSet edges = collectEdgeSet(tree);
    int common = 0;
    for (auto const& e : edges) if (target_edges.count(e)) ++common;
    return common;
}

// Simple shape detector: returns true when every node has at most one child,
// which is characteristic of comb/path-like trees.
bool looksLikeComb(const VectorRangeTreeMap& tree) {
    for (int v : tree.original_nodes) {
        int cnt = 0;
        int L = tree.getLeftChild(v);
        int R = tree.getRightChild(v);
        if (L != VectorRangeTreeMap::NO_CHILD && tree.isOriginal(L)) ++cnt;
        if (R != VectorRangeTreeMap::NO_CHILD && tree.isOriginal(R)) ++cnt;
        if (cnt > 1) return false;
    }
    return true;
}

// Lightweight helper that interprets environment flags such as
// MY_BFS_FILTER_ROTATIONS, returning true when the first character is
// truthy (non-zero, non-f/n).
bool envFlagEnabled(const char* name) {
    if (const char* val = std::getenv(name)) {
        if (val[0] == '\0' || val[0] == '0') return false;
        if (val[0] == 'f' || val[0] == 'F') return false;
        if (val[0] == 'n' || val[0] == 'N') return false;
        return true;
    }
    return false;
}

size_t envSize(const char* name, size_t fallback) {
    if (const char* val = std::getenv(name)) {
        char* end = nullptr;
        errno = 0;
        unsigned long long parsed = std::strtoull(val, &end, 10);
        if (errno == 0 && end && end != val)
            return static_cast<size_t>(parsed);
    }
    return fallback;
}

bool visitedDebugEnabled() { static const bool enabled = envFlagEnabled("BFS_VISITED_DEBUG"); return enabled; }
bool statsPrintEnabled()   { static const bool enabled = envFlagEnabled("BFS_PRINT_STATS");    return enabled; }
bool rotationFilterEnabled() { static const bool enabled = envFlagEnabled("MY_BFS_FILTER_ROTATIONS"); return enabled; }
bool optimizedModeEnabled() {
    static int mode = [](){
        if (const char* val = std::getenv("MY_BFS_MODE")) {
            if (val[0] == '\0') return 1;
            char c = static_cast<char>(std::tolower(static_cast<unsigned char>(val[0])));
            if (c == 'b' || c == '0' || c == 'f' || c == 'n') return 0;
        }
        return 1;
    }();
    return mode != 0;
}
bool bidirFallbackEnabled() {
    static int mode = [](){
        if (const char* val = std::getenv("MY_BFS_USE_BIDIR")) {
            if (val[0] == '\0' || val[0] == '0') return 0;
            char c = static_cast<char>(std::tolower(static_cast<unsigned char>(val[0])));
            if (c == 'f' || c == 'n') return 0;
            if (c == '1' || c == 'y' || c == 't') return 1;
        }
        return 1; // default: enable fallback
    }();
    return mode != 0;
}
bool bestFirstEnabled() {
    static int mode = [](){
        if (const char* val = std::getenv("MY_BFS_HEURISTIC")) {
            if (val[0] == '\0') return 0;
            char c = static_cast<char>(std::tolower(static_cast<unsigned char>(val[0])));
            if (c == 'b') return 1;
        }
        return 0;
    }();
    return mode != 0;
}
bool symmetryPruneEnabled() { static const bool enabled = envFlagEnabled("MY_BFS_SYMMETRY_PRUNE"); return enabled; }
size_t transpositionCap() {
    static size_t cap = envSize("MY_BFS_TRANSPOSITION_CAP", 0);
    return cap;
}

size_t nodeCount(const VectorRangeTreeMap& tree) {
    return tree.original_nodes.size();
}

size_t autoFilterThreshold() {
    static size_t threshold = envSize("MY_BFS_AUTO_FILTER_THRESHOLD", 10);
    return threshold;
}

bool preferBiBfsDefault(const VectorRangeTreeMap& A, const VectorRangeTreeMap& B) {
    static size_t threshold = envSize("MY_BFS_BIDIR_PREFER_THRESHOLD", 12);
    if (looksLikeComb(A) && looksLikeComb(B)) return false;
    return nodeCount(A) >= threshold;
}

bool autoMonotonicEnabled(const VectorRangeTreeMap& A, const VectorRangeTreeMap& B) {
    if (looksLikeComb(A) && looksLikeComb(B)) return false;
    return nodeCount(A) >= autoFilterThreshold();
}
} // namespace

// Baseline BFS implementation (no rotation filtering or parent-hash pruning).
int BFSSearchBaseline(const VectorRangeTreeMap& T_s, const VectorRangeTreeMap& T_e,
                      BFSStats* stats)
{
    BFSStats dummy{};
    BFSStats* diag = stats ? stats : &dummy;
    *diag = BFSStats{};

    if (TreesEqual(T_s, T_e)) {
        diag->equality_hits = 1;
        diag->enqueued = 1;
        diag->visited = 1;
        diag->max_queue = 1;
        return 0;
    }

    std::queue<std::pair<VectorRangeTreeMap,int>> Q;
    std::unordered_set<uint64_t> visited;
    visited.reserve(1 << 16);

    auto push_state = [&](VectorRangeTreeMap&& tree, int depth) {
        uint64_t key = canonicalHash(tree);
        if (!visited.emplace(key).second) {
            ++diag->duplicates;
            return;
        }
        Q.emplace(std::move(tree), depth);
        ++diag->enqueued;
        if (Q.size() > diag->max_queue) diag->max_queue = Q.size();
    };

    push_state(VectorRangeTreeMap(T_s), 0);
    diag->visited = visited.size();

    while (!Q.empty()) {
        auto cur_pair = std::move(Q.front()); Q.pop();
        VectorRangeTreeMap cur = std::move(cur_pair.first);
        int d = cur_pair.second;

        ++diag->generated;
        diag->visited = visited.size();
        const int nd = d + 1;

        for (int v : cur.original_nodes) {
            if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur;
                tmp.rotateLeft(v);
                if (TreesEqual(tmp, T_e)) {
                    ++diag->equality_hits;
                    diag->visited = visited.size();
                    return nd;
                }
                push_state(std::move(tmp), nd);
            }
            if (cur.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                auto tmp = cur;
                tmp.rotateRight(v);
                if (TreesEqual(tmp, T_e)) {
                    ++diag->equality_hits;
                    diag->visited = visited.size();
                    return nd;
                }
                push_state(std::move(tmp), nd);
            }
        }
    }

    diag->visited = visited.size();
    return INT_MAX;
}

// Breadth-first search from T_s to T_e with optional stats collection. The core
// loop hashes each state once, skips obvious undo moves, and (optionally) avoids
// exploring rotations that would reduce the number of target edges (when both
// trees look like combs), with a fallback to the unfiltered expansion when no
// non-decreasing move exists.
int BFSSearchOptimized(const VectorRangeTreeMap& T_s, const VectorRangeTreeMap& T_e,
                       BFSStats* stats)
{
    BFSStats local_stats{};
    BFSStats* diag = stats ? stats : &local_stats;
    *diag = BFSStats{};

    const bool debugVisited = visitedDebugEnabled();
    const bool printSummary = statsPrintEnabled();

    if (TreesEqual(T_s, T_e)) {
        diag->equality_hits = 1;
        diag->enqueued = 1;
        diag->visited = 1;
        diag->max_queue = 1;
        if (printSummary) {
            std::cerr << "[BFSSearch] start == target; enqueued=1 visited=1\n";
        }
        return 0;
    }

    const uint64_t tgt_key = canonicalHash(T_e);

    std::unordered_set<uint64_t> visited;
    visited.reserve(1 << 16);

    std::unordered_map<uint64_t, std::string> debug_signatures;
    if (debugVisited) debug_signatures.reserve(1 << 16);

    auto finish = [&](int result) -> int {
        diag->visited = visited.size();
        if (printSummary) {
            std::cerr << "[BFSSearch] dist=" << result
                      << " generated=" << diag->generated
                      << " enqueued=" << diag->enqueued
                      << " duplicates=" << diag->duplicates
                      << " equality_hits=" << diag->equality_hits
                      << " visited=" << diag->visited
                      << " max_queue=" << diag->max_queue
                      << " mismatches=" << diag->signature_mismatches
                      << "\n";
        }
        return result;
    };

    struct Node {
        VectorRangeTreeMap tree;
        int depth;
        uint64_t hash;
        uint64_t parent_hash;
    };

    const bool use_best_first = bestFirstEnabled();

    std::queue<Node> Q;
    struct QueueEntry {
        int f;
        uint64_t order;
        Node node;
    };
    struct CompareEntry {
        bool operator()(const QueueEntry& a, const QueueEntry& b) const {
            if (a.f != b.f) return a.f > b.f;
            return a.order > b.order;
        }
    };
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, CompareEntry> PQ;
    uint64_t push_counter = 0;

    const size_t tt_cap = transpositionCap();
    std::vector<std::optional<int>> transposition_table;
    if (tt_cap) transposition_table.resize(tt_cap);

    const bool comb_shape = looksLikeComb(T_s) && looksLikeComb(T_e);
    const bool comb_mode = comb_shape;
    const bool want_filter = rotationFilterEnabled();
    const bool auto_monotonic = autoMonotonicEnabled(T_s, T_e);
    const bool monotonic_mode = (want_filter && !comb_shape) || auto_monotonic;
    const bool filter_enabled = comb_mode || monotonic_mode;
    const bool symmetry_enabled = symmetryPruneEnabled();

    EdgeSet target_edges = (filter_enabled || use_best_first) ? collectEdgeSet(T_e) : EdgeSet{};
    const int target_edge_total = static_cast<int>(target_edges.size());

    auto heuristic_for = [&](const VectorRangeTreeMap& tree)->int {
        if (!use_best_first || target_edge_total == 0) return 0;
        int missing = target_edge_total - countMatchingTargetEdges(tree, target_edges);
        return missing < 0 ? 0 : missing;
    };

    auto push_to_queue = [&](Node&& node) {
        int h = heuristic_for(node.tree);
        if (use_best_first) {
            PQ.push(QueueEntry{node.depth + h, push_counter++, std::move(node)});
            diag->max_queue = std::max<uint64_t>(diag->max_queue, PQ.size());
        } else {
            Q.push(std::move(node));
            diag->max_queue = std::max<uint64_t>(diag->max_queue, Q.size());
        }
        ++diag->enqueued;
    };

    auto queue_empty = [&]()->bool {
        return use_best_first ? PQ.empty() : Q.empty();
    };

    auto pop_node = [&]()->Node {
        if (use_best_first) {
            Node node = std::move(PQ.top().node);
            PQ.pop();
            return node;
        }
        Node node = std::move(Q.front());
        Q.pop();
        return node;
    };

    auto enqueue_state = [&](VectorRangeTreeMap&& tmp,
                              int depth,
                              uint64_t parent_hash,
                              const char* tag) -> std::optional<int>
    {
        ++diag->generated;

        uint64_t key = canonicalHash(tmp);

        if (key == parent_hash) {
            return std::nullopt; // immediate undo; skip
        }

        if (key == tgt_key) {
            if (TreesEqual(tmp, T_e)) {
                ++diag->equality_hits;
                return depth;
            }
            ++diag->signature_mismatches;
            if (debugVisited) {
                std::cerr << "[BFSSearch] WARNING: hash matched target but TreesEqual failed via "
                          << tag << "\n";
            }
        }

        auto [_, inserted] = visited.emplace(key);
        if (!inserted) {
            ++diag->duplicates;
            if (debugVisited) {
                std::string sig = canonicalKey(tmp);
                auto it = debug_signatures.find(key);
                if (it == debug_signatures.end()) {
                    debug_signatures.emplace(key, sig);
                    std::cerr << "[BFSSearch] WARNING: duplicate key missing cache via "
                              << tag << "\n";
                } else if (it->second != sig) {
                    ++diag->signature_mismatches;
                    std::cerr << "[BFSSearch] WARNING: hash collision detected via "
                              << tag << "\n";
                }
            }
            return std::nullopt;
        }

        if (debugVisited) {
            debug_signatures.emplace(key, canonicalKey(tmp));
        }

        if (tt_cap) {
            size_t idx = static_cast<size_t>(key % tt_cap);
            if (transposition_table[idx].has_value() && transposition_table[idx].value() <= depth) {
                return std::nullopt;
            }
            transposition_table[idx] = depth;
        }

        push_to_queue(Node{std::move(tmp), depth, key, parent_hash});
        return std::nullopt;
    };

    {
        VectorRangeTreeMap start_copy(T_s);
        uint64_t start_hash = canonicalHash(start_copy);
        visited.emplace(start_hash);
        if (debugVisited) debug_signatures.emplace(start_hash, canonicalKey(start_copy));
        push_to_queue(Node{std::move(start_copy), 0, start_hash,
                    std::numeric_limits<uint64_t>::max()});
    }

    while (!queue_empty()) {
        Node cur = pop_node();
        const int d = cur.depth;
        const int nd = d + 1;

        std::unordered_map<int,uint64_t> subtree_hash_memo;
        auto get_sub_hash = [&](int node)->uint64_t {
            return subtreeHash(cur.tree, node, subtree_hash_memo);
        };

        const int base_matches = monotonic_mode ? countMatchingTargetEdges(cur.tree, target_edges) : 0;
        bool filtered_any = false;
        bool accepted_any = false;

        auto symmetry_skip = [&](int v, bool left)->bool {
            if (!symmetry_enabled) return false;
            int L = cur.tree.getLeftChild(v);
            int R = cur.tree.getRightChild(v);
            if (L == VectorRangeTreeMap::NO_CHILD || R == VectorRangeTreeMap::NO_CHILD) return false;
            if (!cur.tree.isOriginal(L) || !cur.tree.isOriginal(R)) return false;
            return left && get_sub_hash(L) == get_sub_hash(R);
        };

        auto try_rotation = [&](int v, bool left, bool force)->std::optional<int> {
            if (symmetry_skip(v, left)) return std::nullopt;
            if (comb_mode && !force) {
                if (!rotationIntroducesTargetEdge(cur.tree, target_edges, v, left)) {
                    filtered_any = true;
                    return std::nullopt;
                }
            }
            auto tmp = cur.tree;
            if (left) tmp.rotateLeft(v);
            else      tmp.rotateRight(v);

            if (monotonic_mode && !force) {
                int after_matches = countMatchingTargetEdges(tmp, target_edges);
                if (after_matches < base_matches) {
                    filtered_any = true;
                    return std::nullopt;
                }
            }

            accepted_any = true;
            auto res = enqueue_state(std::move(tmp), nd, cur.hash, left ? "rotateLeft" : "rotateRight");
            if (res.has_value()) return res;
            return std::nullopt;
        };

        auto run_pass = [&](bool force)->std::optional<int> {
            for (int v : cur.tree.original_nodes) {
                if (cur.tree.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                    if (auto hit = try_rotation(v, /*left=*/true, force); hit.has_value())
                        return hit;
                }
                if (cur.tree.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                    if (auto hit = try_rotation(v, /*left=*/false, force); hit.has_value())
                        return hit;
                }
            }
            return std::nullopt;
        };

        if (auto hit = run_pass(/*force=*/false); hit.has_value())
            return finish(*hit);
        if (filter_enabled && filtered_any && !accepted_any) {
            accepted_any = false;
            if (auto hit = run_pass(/*force=*/true); hit.has_value())
                return finish(*hit);
        }
    }

    return finish(INT_MAX);
}

int BFSSearch(const VectorRangeTreeMap& T_s, const VectorRangeTreeMap& T_e,
              BFSStats* stats)
{
    if (preferBiBfsDefault(T_s, T_e)) {
        size_t cap = envSize("MY_BFS_BIDIR_CAP", 5'000'000);
        if (stats) *stats = BFSStats{};
        int dist = BiBFSSearchHashed(T_s, T_e, cap);
        if (dist != INT_MAX)
            return dist;
    }
    if (optimizedModeEnabled())
        return BFSSearchOptimized(T_s, T_e, stats);
    return BFSSearchBaseline(T_s, T_e, stats);
}

// Symmetric bidirectional BFS with hashed states and the same pruning/fallback
// rules as BFSSearch. Exposed so callers can opt into a meet-in-the-middle pass.
int BiBFSSearchHashed(const VectorRangeTreeMap& S, const VectorRangeTreeMap& T,
                      size_t state_cap)
{
    if (TreesEqual(S, T)) return 0;

    using TV = VectorRangeTreeMap;
    struct Node {
        TV tree;
        int depth;
        uint64_t hash;
        uint64_t parent_hash;
    };

    std::queue<Node> qS, qT;
    std::unordered_map<uint64_t,int> dS, dT;
    dS.reserve(1<<16); dT.reserve(1<<16);

    const bool debugVisited = visitedDebugEnabled();
    std::unordered_map<uint64_t,std::string> dbgS, dbgT;
    if (debugVisited) { dbgS.reserve(1<<16); dbgT.reserve(1<<16); }

    const bool comb_shape = looksLikeComb(S) && looksLikeComb(T);
    const bool want_filter = rotationFilterEnabled();
    const bool auto_monotonic = autoMonotonicEnabled(S, T);
    const bool comb_mode = comb_shape;
    const bool monotonic_mode = (want_filter && !comb_shape) || auto_monotonic;
    const bool filter_enabled = comb_mode || monotonic_mode;
    const bool symmetry_enabled = symmetryPruneEnabled();
    EdgeSet target_edges_S = filter_enabled ? collectEdgeSet(T) : EdgeSet{}; // from S side, target is T
    EdgeSet target_edges_T = filter_enabled ? collectEdgeSet(S) : EdgeSet{}; // from T side, target is S

    auto insert_root = [&](const VectorRangeTreeMap& root,
                           std::queue<Node>& q,
                           std::unordered_map<uint64_t,int>& dist,
                           std::unordered_map<uint64_t,std::string>& dbg) {
        uint64_t key = canonicalHash(root);
        dist.emplace(key, 0);
        if (debugVisited) dbg.emplace(key, canonicalKey(root));
        q.push(Node{root, 0, key, std::numeric_limits<uint64_t>::max()});
    };

    insert_root(S, qS, dS, dbgS);
    insert_root(T, qT, dT, dbgT);

    auto expand = [&](bool fromS)->std::optional<int> {
        auto& q    = fromS ? qS : qT;
        auto& mine = fromS ? dS : dT;
        auto& othr = fromS ? dT : dS;
        auto& dbgMine = fromS ? dbgS : dbgT;
        auto& dbgOthr = fromS ? dbgT : dbgS;
        EdgeSet& target_edges = fromS ? target_edges_S : target_edges_T;

        size_t layer = q.size();
        while (layer--) {
            Node cur = std::move(q.front());
            q.pop();
            int dcur = cur.depth;

            const int base_matches = (monotonic_mode && filter_enabled)
                                      ? countMatchingTargetEdges(cur.tree, target_edges)
                                      : 0;
            bool filtered_any = false;
            bool accepted_any = false;

            std::unordered_map<int,uint64_t> subtree_hash_memo;
            auto get_sub_hash = [&](int node)->uint64_t {
                return subtreeHash(cur.tree, node, subtree_hash_memo);
            };

            auto symmetry_skip = [&](int v, bool left)->bool {
                if (!symmetry_enabled) return false;
                int L = cur.tree.getLeftChild(v);
                int R = cur.tree.getRightChild(v);
                if (L == TV::NO_CHILD || R == TV::NO_CHILD) return false;
                if (!cur.tree.isOriginal(L) || !cur.tree.isOriginal(R)) return false;
                return left && get_sub_hash(L) == get_sub_hash(R);
            };

            auto consider = [&](int node, bool left_rotation, bool force) -> std::optional<int> {
                    if (symmetry_skip(node, left_rotation)) return std::nullopt;
                    if (comb_mode && filter_enabled && !force) {
                        if (!rotationIntroducesTargetEdge(cur.tree, target_edges, node, left_rotation)) {
                            filtered_any = true;
                            return std::nullopt;
                        }
                    }
                    auto tmp = cur.tree;
                    if (left_rotation) tmp.rotateLeft(node);
                    else               tmp.rotateRight(node);

                    uint64_t key = canonicalHash(tmp);

                    if (key == cur.parent_hash) return std::nullopt;

                    if (monotonic_mode && filter_enabled && !force) {
                        int after_matches = countMatchingTargetEdges(tmp, target_edges);
                        if (after_matches < base_matches) {
                            filtered_any = true;
                            return std::nullopt;
                        }
                    }

                    accepted_any = true;
                    if (auto it = othr.find(key); it != othr.end()) {
                        if (debugVisited) {
                            std::string sig = canonicalKey(tmp);
                            auto jt = dbgOthr.find(key);
                            if (jt == dbgOthr.end()) {
                                dbgOthr.emplace(key, sig);
                                std::cerr << "[BiBFSSearch] WARNING: missing opposite signature cache\n";
                            } else if (jt->second != sig) {
                                std::cerr << "[BiBFSSearch] WARNING: hash collision detected between frontiers\n";
                                return std::nullopt;
                            }
                        }
                        return dcur + 1 + it->second;
                    }

                    auto [_, inserted] = mine.emplace(key, dcur + 1);
                    if (!inserted) return std::nullopt;

                    if (debugVisited) dbgMine.emplace(key, canonicalKey(tmp));

                    q.push(Node{std::move(tmp), dcur + 1, key, cur.hash});
                    return std::nullopt;
                };

            auto run_pass = [&](bool force)->std::optional<int> {
                for (int v : cur.tree.original_nodes) {
                    if (cur.tree.getRightChild(v) != TV::NO_CHILD) {
                        if (auto hit = consider(v, true, force); hit.has_value())
                            return hit;
                    }
                    if (cur.tree.getLeftChild(v) != TV::NO_CHILD) {
                        if (auto hit = consider(v, false, force); hit.has_value())
                            return hit;
                    }
                }
                return std::nullopt;
            };

            if (auto hit = run_pass(/*force=*/false); hit.has_value())
                return hit;
            if (filter_enabled && filtered_any && !accepted_any) {
                accepted_any = false;
                if (auto hit = run_pass(/*force=*/true); hit.has_value())
                    return hit;
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

// A feature-rich BFS wrapper that exposes timing/queue/visited caps and records
// telemetry for benchmarking. Shares the same pruning/fallback logic as BFSSearch.
BFSRun BFSSearchCapped(const VectorRangeTreeMap& T_s,
                       const VectorRangeTreeMap& T_e,
                       double time_limit_sec,
                       size_t visited_cap,
                       size_t queue_cap)
{
    using Clock = std::chrono::steady_clock;
    const auto t0 = Clock::now();

    const uint64_t tgt_key = canonicalHash(T_e);

    struct Node {
        VectorRangeTreeMap tree;
        int depth;
        uint64_t hash;
        uint64_t parent_hash;
    };

    std::queue<Node> Q;
    std::unordered_set<uint64_t> V;
    V.reserve(std::min(visited_cap, size_t(2'000'000)));

    size_t expanded = 0;
    size_t enqueued = 0;
    size_t duplicates = 0;
    size_t equality_hits = 0;
    size_t signature_mismatches = 0;
    size_t max_queue = 0;
    size_t generated = 0;

    const bool debugVisited = visitedDebugEnabled();
    const bool printSummary = statsPrintEnabled();

    std::unordered_map<uint64_t, std::string> debug_hash_map;
    if (debugVisited) debug_hash_map.reserve(1 << 16);

    auto elapsed_seconds = [&]() -> double {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };

    auto make_result = [&](int dist, bool timeout, bool cap_hit, double elapsed) -> BFSRun {
        BFSRun res{dist, expanded, enqueued, V.size(), timeout, cap_hit, elapsed,
                   duplicates, max_queue, equality_hits, signature_mismatches};
        if (printSummary) {
            std::cerr << "[BFSSearchCapped] dist=" << dist
                      << " expanded=" << expanded
                      << " enqueued=" << enqueued
                      << " duplicates=" << duplicates
                      << " equality_hits=" << equality_hits
                      << " max_queue=" << max_queue
                      << " visited=" << res.visited
                      << " mismatches=" << signature_mismatches
                      << " generated=" << generated
                      << (timeout ? " TIMEOUT" : "")
                      << (cap_hit ? " CAP" : "")
                      << " seconds=" << elapsed
                      << "\n";
        }
        return res;
    };

    auto try_bidir = [&]() -> std::optional<BFSRun> {
        if (!bidirFallbackEnabled()) return std::nullopt;
        size_t cap = envSize("MY_BFS_BIDIR_CAP", 5'000'000);
        int d_bidir = BiBFSSearchHashed(T_s, T_e, cap);
        if (d_bidir != INT_MAX)
            return make_result(d_bidir, false, false, elapsed_seconds());
        return std::nullopt;
    };

    const bool want_filter = rotationFilterEnabled();
    const bool comb_mode = want_filter && looksLikeComb(T_s) && looksLikeComb(T_e);
    const bool monotonic_mode = want_filter && !comb_mode;
    const bool filter_enabled = comb_mode || monotonic_mode;
    EdgeSet target_edges = filter_enabled ? collectEdgeSet(T_e) : EdgeSet{};

    auto attempt_enqueue = [&](VectorRangeTreeMap&& X,
                               int depth,
                               uint64_t parent_hash,
                               const char* tag) -> std::optional<BFSRun>
    {
        ++generated;

        uint64_t key = canonicalHash(X);

        if (key == parent_hash) {
            return std::nullopt;
        }

        if (key == tgt_key) {
            if (TreesEqual(X, T_e)) {
                ++equality_hits;
                return make_result(depth, false, false, elapsed_seconds());
            }
            ++signature_mismatches;
            if (debugVisited) {
                std::cerr << "[BFSSearchCapped] WARNING: hash matched target but TreesEqual failed via "
                          << tag << "\n";
            }
        }

        if (V.size() >= visited_cap || Q.size() >= queue_cap) {
            if (auto fb = try_bidir(); fb.has_value())
                return *fb;
            return make_result(INT_MAX, false, true, elapsed_seconds());
        }

        auto [_, inserted] = V.emplace(key);
        if (!inserted) {
            ++duplicates;
            if (debugVisited) {
                std::string sig = canonicalKey(X);
                auto it = debug_hash_map.find(key);
                if (it == debug_hash_map.end()) {
                    debug_hash_map.emplace(key, sig);
                    std::cerr << "[BFSSearchCapped] WARNING: duplicate key missing reference via "
                              << tag << "\n";
                } else if (it->second != sig) {
                    ++signature_mismatches;
                    std::cerr << "[BFSSearchCapped] WARNING: hash collision detected via "
                              << tag << "\n";
                }
            }
            return std::nullopt;
        }

        if (debugVisited) {
            debug_hash_map.emplace(key, canonicalKey(X));
        }

        Q.push(Node{std::move(X), depth, key, parent_hash});
        ++enqueued;
        if (Q.size() > max_queue) max_queue = Q.size();
        return std::nullopt;
    };

    {
        VectorRangeTreeMap start_copy(T_s);
        uint64_t start_hash = canonicalHash(start_copy);
        V.emplace(start_hash);
        if (debugVisited) debug_hash_map.emplace(start_hash, canonicalKey(start_copy));
        Q.push(Node{std::move(start_copy), 0, start_hash,
                    std::numeric_limits<uint64_t>::max()});
        enqueued = 1;
        max_queue = 1;
    }

    while (!Q.empty()) {
        const double sec = std::chrono::duration<double>(Clock::now() - t0).count();
        if (sec > time_limit_sec) {
            if (auto fb = try_bidir(); fb.has_value())
                return *fb;
            return make_result(INT_MAX, true, false, sec);
        }

        Node cur = std::move(Q.front());
        Q.pop();
        int d     = cur.depth;

        ++expanded;
        const int nextD = d + 1;

        const int base_matches = monotonic_mode ? countMatchingTargetEdges(cur.tree, target_edges) : 0;
        bool filtered_any = false;
        bool accepted_any = false;

        auto try_rotation = [&](int v, bool left, bool force)->std::optional<BFSRun> {
            if (comb_mode && !force) {
                if (!rotationIntroducesTargetEdge(cur.tree, target_edges, v, left)) {
                    filtered_any = true;
                    return std::nullopt;
                }
            }
            auto tmp = cur.tree;
            if (left) tmp.rotateLeft(v);
            else      tmp.rotateRight(v);

            if (monotonic_mode && !force) {
                int after_matches = countMatchingTargetEdges(tmp, target_edges);
                if (after_matches < base_matches) {
                    filtered_any = true;
                    return std::nullopt;
                }
            }

            accepted_any = true;
            auto res = attempt_enqueue(std::move(tmp), nextD, cur.hash, left ? "rotateLeft" : "rotateRight");
            if (res.has_value()) return res;
            return std::nullopt;
        };

        auto run_pass = [&](bool force)->std::optional<BFSRun> {
            for (int v : cur.tree.original_nodes) {
                if (cur.tree.getRightChild(v) != VectorRangeTreeMap::NO_CHILD) {
                    if (auto res = try_rotation(v, /*left=*/true, force); res.has_value())
                        return res;
                }
                if (cur.tree.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD) {
                    if (auto res = try_rotation(v, /*left=*/false, force); res.has_value())
                        return res;
                }
            }
            return std::nullopt;
        };

        if (auto res = run_pass(/*force=*/false); res.has_value())
            return *res;
        if (filter_enabled && filtered_any && !accepted_any) {
            accepted_any = false;
            if (auto res = run_pass(/*force=*/true); res.has_value())
                return *res;
        }
    }

    if (auto fb = try_bidir(); fb.has_value())
        return *fb;
    return make_result(INT_MAX, false, false, elapsed_seconds());
}
