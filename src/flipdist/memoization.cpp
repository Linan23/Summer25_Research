#include "memoization.h"
#include "algorithm.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// Debug flag
const bool DEBUG = false;

ProfileStats g_profile;

std::unordered_map<Key128, bool, Key128Hash> g_memo_flipdist;
std::unordered_map<Key128, bool, Key128Hash> g_memo_tdi;
std::unordered_map<Key128, bool, Key128Hash> g_memo_tds;
std::unordered_map<Key128, KBounds, Key128Hash> g_kbounds;
std::unordered_map<Key128, KBounds, Key128Hash> g_tds_bounds;
std::unordered_map<Key128, KBounds, Key128Hash> g_tdi_bounds;
static std::unordered_map<Key128, PlateauStateRecord, Key128Hash> g_plateau_states;
static bool g_plateau_profile_enabled = false;
static std::string g_plateau_output_path;
static std::string g_plateau_case_type;
static std::string g_plateau_direction;
static int g_plateau_n = -1;
static int g_plateau_seed = -1;
static bool g_plateau_output_dirty = false;
static std::uint64_t g_plateau_snapshot_seq = 0;
static int g_plateau_samples_since_flush = 0;
static std::chrono::steady_clock::time_point g_plateau_last_flush;
static std::unordered_map<Key128, TdiPostIStateRecord, Key128Hash> g_tdi_posti_states;
static bool g_tdi_posti_profile_enabled = false;
static std::string g_tdi_posti_output_path;
static std::string g_tdi_posti_case_type;
static std::string g_tdi_posti_direction;
static int g_tdi_posti_n = -1;
static int g_tdi_posti_seed = -1;

static int plateauFlushEverySamples() {
    const char *env = std::getenv("FLIPDIST_PROFILE_PLATEAU_FLUSH_EVERY");
    if (!env) return 256;
    int value = std::atoi(env);
    return value > 0 ? value : 256;
}

static int plateauFlushIntervalMs() {
    const char *env = std::getenv("FLIPDIST_PROFILE_PLATEAU_FLUSH_MS");
    if (!env) return 500;
    int value = std::atoi(env);
    return value > 0 ? value : 500;
}

size_t Key128Hash::operator()(const Key128& k) const {
    uint64_t h = k.hi;
    h ^= (k.lo + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
    return static_cast<size_t>(h);
}

struct TreeFingerprint {
    uint64_t h1 = 0;
    uint64_t h2 = 0;
};

static uint64_t fnvMix(uint64_t h, uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
    return h;
}

static uint64_t edgeKey(const std::pair<int,int>& e) {
    int a = std::min(e.first, e.second);
    int b = std::max(e.first, e.second);
    return (static_cast<uint64_t>(a) << 32) ^ static_cast<uint32_t>(b);
}

static std::pair<int,int> canonicalEdge(std::pair<int,int> e) {
    if (e.second < e.first) {
        std::swap(e.first, e.second);
    }
    return e;
}

static std::pair<std::pair<int,int>, std::pair<int,int>> canonicalPair(
    std::pair<std::pair<int,int>, std::pair<int,int>> p) {
    p.first = canonicalEdge(p.first);
    p.second = canonicalEdge(p.second);
    if (p.second < p.first) {
        std::swap(p.first, p.second);
    }
    return p;
}

static uint64_t hashEdgeSet(const std::vector<std::pair<int,int>>& edges) {
    std::vector<uint64_t> keys;
    keys.reserve(edges.size());
    for (const auto& e : edges) keys.push_back(edgeKey(e));
    std::sort(keys.begin(), keys.end());
    uint64_t h = 1469598103934665603ULL;
    for (uint64_t v : keys) {
        h = fnvMix(h, v);
    }
    return h;
}

static uint64_t hashPairSet(const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& pairs) {
    std::vector<std::pair<uint64_t,uint64_t>> keys;
    keys.reserve(pairs.size());
    for (const auto& p : pairs) {
        uint64_t k1 = edgeKey(p.first);
        uint64_t k2 = edgeKey(p.second);
        if (k2 < k1) std::swap(k1, k2);
        keys.emplace_back(k1, k2);
    }
    std::sort(keys.begin(), keys.end());
    uint64_t h = 1469598103934665603ULL;
    for (const auto& kv : keys) {
        h = fnvMix(h, kv.first);
        h = fnvMix(h, kv.second);
    }
    return h;
}

void dedupeSPairs(std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& pairs) {
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> out;
    out.reserve(pairs.size());
    std::unordered_set<Key128, Key128Hash> seen;
    for (const auto& p : pairs) {
        auto cp = canonicalPair(p);
        uint64_t k1 = edgeKey(cp.first);
        uint64_t k2 = edgeKey(cp.second);
        if (k2 < k1) std::swap(k1, k2);
        Key128 key{k1, k2};
        if (seen.insert(key).second) {
            out.push_back(cp);
        }
    }
    pairs.swap(out);
}

static void mixKey(Key128& k, uint64_t v) {
    k.hi = fnvMix(k.hi, v);
    k.lo = fnvMix(k.lo, v ^ 0x9e3779b97f4a7c15ULL);
}

static void mixFingerprint(TreeFingerprint& fp, uint64_t v) {
    fp.h1 = fnvMix(fp.h1, v);
    fp.h2 = fnvMix(fp.h2, v ^ 0x9e3779b97f4a7c15ULL);
}

static TreeFingerprint treeFingerprint(const VectorRangeTreeMap& T) {
    if (T.fingerprint_valid) {
        return TreeFingerprint{T.fingerprint_h1, T.fingerprint_h2};
    }
    TreeFingerprint out{1469598103934665603ULL, 1099511628211ULL};
    mixFingerprint(out, static_cast<uint64_t>(T.original_inorder.size()));
    mixFingerprint(out, static_cast<uint64_t>(T.root + 1));
    for (int v : T.original_inorder) {
        const bool in_range = (v >= 0 && v < static_cast<int>(T.ranges.size()) &&
                               v < static_cast<int>(T.edges.size()) &&
                               v < static_cast<int>(T.parents.size()));
        const auto r = in_range ? T.ranges[v] : std::pair<int, int>{0, 0};
        const int L = in_range ? T.edges[v].first : VectorRangeTreeMap::NO_CHILD;
        const int R = in_range ? T.edges[v].second : VectorRangeTreeMap::NO_CHILD;
        const int P = in_range ? T.parents[v] : VectorRangeTreeMap::NO_PARENT;
        mixFingerprint(out, static_cast<uint64_t>(v + 1));
        mixFingerprint(out, (static_cast<uint64_t>(static_cast<uint32_t>(r.first)) << 32) ^
                            static_cast<uint32_t>(r.second));
        mixFingerprint(out, (static_cast<uint64_t>(static_cast<uint32_t>(L + 2)) << 32) ^
                            static_cast<uint32_t>(R + 2));
        mixFingerprint(out, static_cast<uint64_t>(P + 2));
    }
    T.fingerprint_h1 = out.h1;
    T.fingerprint_h2 = out.h2;
    T.fingerprint_valid = true;
    return out;
}

static Key128 makeTreePairKeyFromFingerprints(int tag, const TreeFingerprint& a,
                                              const TreeFingerprint& b) {
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, static_cast<uint64_t>(tag));
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    return key;
}

Key128 makeKeyPair(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    return makeTreePairKeyFromFingerprints(1, a, b);
}

Key128 makeKeyFlip(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 2);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, static_cast<uint64_t>(k));
    return key;
}

Key128 makeKeyI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
                       const std::vector<std::pair<int,int>>& I) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    uint64_t h = hashEdgeSet(I);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 3);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(I.size()));
    mixKey(key, static_cast<uint64_t>(k));
    return key;
}

Key128 makeKeyS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
                       const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    uint64_t h = hashPairSet(S);
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_end);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 4);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(S.size()));
    mixKey(key, static_cast<uint64_t>(k));
    return key;
}

Key128 makeKeySBase(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end,
                           const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    uint64_t h = hashPairSet(S);
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_end);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 5);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(S.size()));
    return key;
}

EmptySKeys makeEmptySKeys(const VectorRangeTreeMap& T_init,
                          const VectorRangeTreeMap& T_end,
                          int k) {
    const TreeFingerprint a = treeFingerprint(T_init);
    const TreeFingerprint b = treeFingerprint(T_end);
    constexpr uint64_t empty_s_hash = 1469598103934665603ULL;
    EmptySKeys keys;
    keys.pair_key = makeTreePairKeyFromFingerprints(1, a, b);

    keys.bounds_key = makeTreePairKeyFromFingerprints(5, a, b);
    mixKey(keys.bounds_key, empty_s_hash);
    mixKey(keys.bounds_key, 0);

    keys.exact_key = makeTreePairKeyFromFingerprints(4, a, b);
    mixKey(keys.exact_key, empty_s_hash);
    mixKey(keys.exact_key, 0);
    mixKey(keys.exact_key, static_cast<uint64_t>(k));
    return keys;
}

Key128 makeKeyIBase(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final,
                           const std::vector<std::pair<int,int>>& I) {
    uint64_t h = hashEdgeSet(I);
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 6);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    mixKey(key, h);
    mixKey(key, static_cast<uint64_t>(I.size()));
    return key;
}

bool tryBoundsPrune(const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                           const Key128& key, int k, bool& value_out) {
    auto it = bounds_map.find(key);
    if (it == bounds_map.end()) return false;
    const auto& b = it->second;
    if (b.min_success >= 0 && k >= b.min_success) {
        value_out = true;
        return true;
    }
    if (k <= b.max_fail) {
        value_out = false;
        return true;
    }
    return false;
}

void updateBounds(std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                         const Key128& key, int k, bool value) {
    auto& b = bounds_map[key];
    bool changed = false;
    if (value) {
        if (b.min_success < 0 || k < b.min_success) {
            b.min_success = k;
            changed = true;
        }
        const int forced_fail = k - 1;
        if (forced_fail < b.max_fail) {
            b.max_fail = forced_fail;
            changed = true;
        }
        if (b.max_fail >= b.min_success) {
            const int needed_max_fail = b.min_success - 1;
            if (b.max_fail != needed_max_fail) {
                b.max_fail = needed_max_fail;
                changed = true;
            }
        }
    } else {
        if (k > b.max_fail) {
            b.max_fail = k;
            changed = true;
        }
    }
    (void)bounds_map;
    (void)key;
    (void)k;
    (void)value;
}

int requiredBudgetFromBounds(const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                    const Key128& key) {
    auto it = bounds_map.find(key);
    if (it == bounds_map.end()) return 0;
    return std::max(0, it->second.max_fail + 1);
}

void resetMemo() {
    g_memo_flipdist.clear();
    g_memo_tdi.clear();
    g_memo_tds.clear();
    g_kbounds.clear();
    g_tds_bounds.clear();
    g_tdi_bounds.clear();
    resetAlgorithmCaches();
}

void initProfile() {
    if (g_profile.enabled) return;
    const char *env = std::getenv("FLIPDIST_PROFILE");
    if (env && std::string(env) == "1") {
        g_profile.enabled = true;
        g_profile.start_time = std::chrono::steady_clock::now();
        const char *plateau_env = std::getenv("FLIPDIST_PROFILE_PLATEAU");
        g_plateau_profile_enabled = (plateau_env && std::string(plateau_env) == "1");
        const char *tdi_posti_env = std::getenv("FLIPDIST_PROFILE_TDI_POSTI");
        g_tdi_posti_profile_enabled = (tdi_posti_env && std::string(tdi_posti_env) == "1");
        const char *abort_env = std::getenv("FLIPDIST_PROFILE_ABORT_MS");
        if (abort_env) {
            g_profile.abort_ms = std::atoi(abort_env);
        }
    }
}

bool plateauProfileEnabled() {
    return g_profile.enabled && g_plateau_profile_enabled;
}

bool tdiPostIProfileEnabled() {
    return g_profile.enabled && g_tdi_posti_profile_enabled;
}

void resetPlateauProfile() {
    g_plateau_states.clear();
    g_plateau_output_dirty = false;
    g_plateau_snapshot_seq = 0;
    g_plateau_samples_since_flush = 0;
    g_plateau_case_type.clear();
    g_plateau_direction.clear();
    g_plateau_n = -1;
    g_plateau_seed = -1;
    g_plateau_last_flush = std::chrono::steady_clock::now();
}

void resetTdiPostIProfile() {
    g_tdi_posti_states.clear();
    g_tdi_posti_case_type.clear();
    g_tdi_posti_direction.clear();
    g_tdi_posti_n = -1;
    g_tdi_posti_seed = -1;
}

static void writePlateauProfileJsonl(std::ostream &os, const std::string &case_type, int n, int seed,
                                     const std::string &direction, const std::string &record_type,
                                     std::uint64_t snapshot_seq) {
    if (!plateauProfileEnabled()) return;
    std::vector<const PlateauStateRecord*> records;
    records.reserve(g_plateau_states.size());
    for (const auto &kv : g_plateau_states) {
        records.push_back(&kv.second);
    }
    std::sort(records.begin(), records.end(), [](const PlateauStateRecord *a, const PlateauStateRecord *b) {
        if (a->fail_count != b->fail_count) return a->fail_count > b->fail_count;
        if (a->recurrence_count != b->recurrence_count) return a->recurrence_count > b->recurrence_count;
        if (a->internal_edges != b->internal_edges) return a->internal_edges > b->internal_edges;
        if (a->state_key.hi != b->state_key.hi) return a->state_key.hi < b->state_key.hi;
        return a->state_key.lo < b->state_key.lo;
    });
    for (const PlateauStateRecord *rec : records) {
        os << std::fixed << std::setprecision(0);
        os << "{"
           << "\"record_type\":\"" << record_type << "\""
           << ",\"snapshot_seq\":" << snapshot_seq
           << ",\"case_type\":\"" << case_type << "\""
           << ",\"n\":" << n
           << ",\"seed\":" << seed
           << ",\"direction\":\"" << direction << "\""
           << ",\"state_hi\":\"" << rec->state_key.hi << "\""
           << ",\"state_lo\":\"" << rec->state_key.lo << "\""
           << ",\"internal_edges\":" << rec->internal_edges
           << ",\"reduced_core_internal_edges\":" << rec->reduced_core_internal_edges
           << ",\"conflicts\":" << rec->conflicts
           << ",\"min_k\":" << rec->min_k
           << ",\"max_k\":" << rec->max_k
           << ",\"legal_children\":" << rec->legal_children
           << ",\"plateau_buckets\":" << rec->plateau_buckets
           << ",\"start_branching_nodes\":" << rec->start_branching_nodes
           << ",\"target_branching_nodes\":" << rec->target_branching_nodes
           << ",\"start_max_branch_subtree_edges\":" << rec->start_max_branch_subtree_edges
           << ",\"target_max_branch_subtree_edges\":" << rec->target_max_branch_subtree_edges
           << ",\"chain_arm_count\":" << rec->chain_arm_count
           << ",\"broom_arm_count\":" << rec->broom_arm_count
           << ",\"branchy_reduction_hits\":" << rec->branchy_reduction_hits
           << ",\"branchy_forced_cost_max\":" << rec->branchy_forced_cost_max
           << ",\"branchy_residual_edges_min\":" << rec->branchy_residual_edges_min
           << ",\"branchy_exact_bound_hits\":" << rec->branchy_exact_bound_hits
           << ",\"branchy_exact_decision_hits\":" << rec->branchy_exact_decision_hits
           << ",\"recurrence_count\":" << rec->recurrence_count
           << ",\"success_count\":" << rec->success_count
           << ",\"fail_count\":" << rec->fail_count
           << ",\"timeout_count\":" << rec->timeout_count
           << ",\"free_edge_hits\":" << rec->free_edge_hits
           << ",\"free_edge_misses\":" << rec->free_edge_misses
           << ",\"start_degree_histogram\":\"" << rec->start_degree_histogram << "\""
           << ",\"target_degree_histogram\":\"" << rec->target_degree_histogram << "\""
           << ",\"motif\":\"" << rec->motif << "\""
           << ",\"core_signature\":\"" << rec->core_signature << "\""
           << ",\"branchy_residual_motif\":\"" << rec->branchy_residual_motif << "\""
           << "}\n";
    }
}

static void flushPlateauProfileToFile(bool final_flush) {
    if (!plateauProfileEnabled() || g_plateau_output_path.empty() || g_plateau_states.empty()) return;
    if (!final_flush && !g_plateau_output_dirty) return;
    std::ofstream out(g_plateau_output_path, std::ios::app);
    if (!out) return;
    const std::string record_type = final_flush ? "plateau_state" : "plateau_state_snapshot";
    const std::uint64_t snapshot_seq = ++g_plateau_snapshot_seq;
    writePlateauProfileJsonl(out, g_plateau_case_type, g_plateau_n, g_plateau_seed,
                             g_plateau_direction, record_type, snapshot_seq);
    out.flush();
    if (!final_flush) {
        g_plateau_output_dirty = false;
        g_plateau_samples_since_flush = 0;
        g_plateau_last_flush = std::chrono::steady_clock::now();
    }
}

static void maybeFlushPlateauProfileSnapshot() {
    if (!plateauProfileEnabled() || g_plateau_output_path.empty()) return;
    if (!g_plateau_output_dirty) return;
    const auto now = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(now - g_plateau_last_flush).count();
    if (g_plateau_samples_since_flush >= plateauFlushEverySamples() ||
        elapsed_ms >= plateauFlushIntervalMs()) {
        flushPlateauProfileToFile(false);
    }
}

void beginPlateauProfileRun(const std::string &case_type, int n, int seed,
                            const std::string &direction, bool clear_output_file) {
    if (!plateauProfileEnabled()) return;
    const char *path_env = std::getenv("FLIPDIST_PROFILE_PLATEAU_FILE");
    g_plateau_output_path = path_env ? path_env : "";
    g_plateau_case_type = case_type;
    g_plateau_n = n;
    g_plateau_seed = seed;
    g_plateau_direction = direction;
    g_plateau_output_dirty = false;
    g_plateau_snapshot_seq = 0;
    g_plateau_samples_since_flush = 0;
    g_plateau_last_flush = std::chrono::steady_clock::now();
    if (clear_output_file && !g_plateau_output_path.empty()) {
        std::ofstream clear(g_plateau_output_path, std::ios::trunc);
    }
}

void beginTdiPostIProfileRun(const std::string &case_type, int n, int seed,
                             const std::string &direction, bool clear_output_file) {
    if (!tdiPostIProfileEnabled()) return;
    const char *path_env = std::getenv("FLIPDIST_PROFILE_TDI_POSTI_FILE");
    g_tdi_posti_output_path = path_env ? path_env : "";
    g_tdi_posti_case_type = case_type;
    g_tdi_posti_n = n;
    g_tdi_posti_seed = seed;
    g_tdi_posti_direction = direction;
    if (clear_output_file && !g_tdi_posti_output_path.empty()) {
        std::ofstream clear(g_tdi_posti_output_path, std::ios::trunc);
    }
}

void notePlateauStateSample(const Key128 &state_key, int internal_edges, int reduced_core_internal_edges,
                            int k, int conflicts, int legal_children, int plateau_buckets,
                            bool free_edge_hit, int start_branching_nodes,
                            int target_branching_nodes, int start_max_branch_subtree_edges,
                            int target_max_branch_subtree_edges, int chain_arm_count,
                            int broom_arm_count, const std::string &start_degree_histogram,
                            const std::string &target_degree_histogram,
                            const std::string &motif, const std::string &core_signature,
                            bool branchy_reduction_triggered, int branchy_forced_cost,
                            int branchy_residual_edges, const std::string &branchy_residual_motif,
                            bool branchy_exact_bound_used, bool branchy_exact_decision_used) {
    if (!plateauProfileEnabled()) return;
    auto &rec = g_plateau_states[state_key];
    rec.state_key = state_key;
    rec.internal_edges = internal_edges;
    rec.reduced_core_internal_edges = reduced_core_internal_edges;
    rec.conflicts = conflicts;
    rec.legal_children = legal_children;
    rec.plateau_buckets = plateau_buckets;
    rec.start_branching_nodes = start_branching_nodes;
    rec.target_branching_nodes = target_branching_nodes;
    rec.start_max_branch_subtree_edges = start_max_branch_subtree_edges;
    rec.target_max_branch_subtree_edges = target_max_branch_subtree_edges;
    rec.chain_arm_count = chain_arm_count;
    rec.broom_arm_count = broom_arm_count;
    if (branchy_reduction_triggered) {
        rec.branchy_reduction_hits++;
        rec.branchy_forced_cost_max = std::max(rec.branchy_forced_cost_max, branchy_forced_cost);
        if (branchy_residual_edges >= 0) {
            if (rec.branchy_residual_edges_min < 0) {
                rec.branchy_residual_edges_min = branchy_residual_edges;
            } else {
                rec.branchy_residual_edges_min =
                    std::min(rec.branchy_residual_edges_min, branchy_residual_edges);
            }
        }
        if (!branchy_residual_motif.empty()) {
            rec.branchy_residual_motif = branchy_residual_motif;
        }
    }
    if (branchy_exact_bound_used) rec.branchy_exact_bound_hits++;
    if (branchy_exact_decision_used) rec.branchy_exact_decision_hits++;
    rec.recurrence_count++;
    if (rec.min_k < 0 || k < rec.min_k) rec.min_k = k;
    if (rec.max_k < 0 || k > rec.max_k) rec.max_k = k;
    if (free_edge_hit) rec.free_edge_hits++;
    else rec.free_edge_misses++;
    if (!motif.empty()) rec.motif = motif;
    if (!core_signature.empty()) rec.core_signature = core_signature;
    if (!start_degree_histogram.empty()) rec.start_degree_histogram = start_degree_histogram;
    if (!target_degree_histogram.empty()) rec.target_degree_histogram = target_degree_histogram;
    g_plateau_output_dirty = true;
    g_plateau_samples_since_flush++;
    maybeFlushPlateauProfileSnapshot();
}

void notePlateauStateOutcome(const Key128 &state_key, const std::string &outcome) {
    if (!plateauProfileEnabled()) return;
    auto it = g_plateau_states.find(state_key);
    if (it == g_plateau_states.end()) return;
    if (outcome == "success") it->second.success_count++;
    else if (outcome == "timeout") it->second.timeout_count++;
    else it->second.fail_count++;
    g_plateau_output_dirty = true;
    maybeFlushPlateauProfileSnapshot();
}

void noteTdiPostIPruneSample(const Key128 &state_key, int remaining_k, int lower_bound,
                             int conflicts, int s_size, int start_branching_nodes,
                             int target_branching_nodes,
                             const std::string &start_degree_histogram,
                             const std::string &target_degree_histogram,
                             const std::string &start_tree,
                             const std::string &target_tree) {
    if (!tdiPostIProfileEnabled()) return;
    auto &rec = g_tdi_posti_states[state_key];
    rec.state_key = state_key;
    rec.conflicts = conflicts;
    rec.s_size = s_size;
    rec.start_branching_nodes = start_branching_nodes;
    rec.target_branching_nodes = target_branching_nodes;
    rec.start_degree_histogram = start_degree_histogram;
    rec.target_degree_histogram = target_degree_histogram;
    rec.start_tree = start_tree;
    rec.target_tree = target_tree;
    rec.recurrence_count++;
    if (rec.min_remaining_k < 0 || remaining_k < rec.min_remaining_k) rec.min_remaining_k = remaining_k;
    if (rec.max_remaining_k < 0 || remaining_k > rec.max_remaining_k) rec.max_remaining_k = remaining_k;
    if (rec.min_lower_bound < 0 || lower_bound < rec.min_lower_bound) rec.min_lower_bound = lower_bound;
    if (rec.max_lower_bound < 0 || lower_bound > rec.max_lower_bound) rec.max_lower_bound = lower_bound;
}

void emitPlateauProfileJsonl(std::ostream &os, const std::string &case_type, int n, int seed,
                             const std::string &direction) {
    writePlateauProfileJsonl(os, case_type, n, seed, direction, "plateau_state", 0);
}

void finalizePlateauProfileRun() {
    if (!plateauProfileEnabled()) return;
    flushPlateauProfileToFile(true);
    g_plateau_output_dirty = false;
}

void emitTdiPostIProfileJsonl(std::ostream &os, const std::string &case_type, int n, int seed,
                              const std::string &direction) {
    if (!tdiPostIProfileEnabled()) return;
    std::vector<const TdiPostIStateRecord*> records;
    records.reserve(g_tdi_posti_states.size());
    for (const auto &kv : g_tdi_posti_states) {
        records.push_back(&kv.second);
    }
    std::sort(records.begin(), records.end(), [](const TdiPostIStateRecord *a, const TdiPostIStateRecord *b) {
        if (a->recurrence_count != b->recurrence_count) return a->recurrence_count > b->recurrence_count;
        if (a->s_size != b->s_size) return a->s_size > b->s_size;
        if (a->state_key.hi != b->state_key.hi) return a->state_key.hi < b->state_key.hi;
        return a->state_key.lo < b->state_key.lo;
    });
    for (const TdiPostIStateRecord *rec : records) {
        os << std::fixed << std::setprecision(0);
        os << "{"
           << "\"record_type\":\"tdi_post_i_state\""
           << ",\"case_type\":\"" << case_type << "\""
           << ",\"n\":" << n
           << ",\"seed\":" << seed
           << ",\"direction\":\"" << direction << "\""
           << ",\"state_hi\":\"" << rec->state_key.hi << "\""
           << ",\"state_lo\":\"" << rec->state_key.lo << "\""
           << ",\"conflicts\":" << rec->conflicts
           << ",\"s_size\":" << rec->s_size
           << ",\"min_remaining_k\":" << rec->min_remaining_k
           << ",\"max_remaining_k\":" << rec->max_remaining_k
           << ",\"min_lower_bound\":" << rec->min_lower_bound
           << ",\"max_lower_bound\":" << rec->max_lower_bound
           << ",\"start_branching_nodes\":" << rec->start_branching_nodes
           << ",\"target_branching_nodes\":" << rec->target_branching_nodes
           << ",\"recurrence_count\":" << rec->recurrence_count
           << ",\"start_degree_histogram\":\"" << rec->start_degree_histogram << "\""
           << ",\"target_degree_histogram\":\"" << rec->target_degree_histogram << "\""
           << ",\"start_tree\":\"" << rec->start_tree << "\""
           << ",\"target_tree\":\"" << rec->target_tree << "\""
           << "}\n";
    }
}

void finalizeTdiPostIProfileRun() {
    if (!tdiPostIProfileEnabled()) return;
    if (!g_tdi_posti_output_path.empty() && !g_tdi_posti_states.empty()) {
        std::ofstream out(g_tdi_posti_output_path, std::ios::app);
        if (out) {
            emitTdiPostIProfileJsonl(out, g_tdi_posti_case_type, g_tdi_posti_n,
                                     g_tdi_posti_seed, g_tdi_posti_direction);
            out.flush();
        }
    }
    resetTdiPostIProfile();
}

bool profileAbortRequested() {
    if (g_profile.abort_ms <= 0) return false;
    auto now = std::chrono::steady_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(now - g_profile.start_time).count();
    return elapsed_ms > g_profile.abort_ms;
}

void throwIfProfileAbort() {
    if (profileAbortRequested()) {
        throw SearchAbortException{};
    }
}

ScopedTimer::ScopedTimer(double *a) : acc(a) {
    if (g_profile.enabled && acc) {
        t0 = std::chrono::steady_clock::now();
    }
}

ScopedTimer::~ScopedTimer() {
    if (g_profile.enabled && acc) {
        auto t1 = std::chrono::steady_clock::now();
        *acc += std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
}

// Definition: Emit a debug line when DEBUG is enabled
// Parameters: msg: message to print
// Returns: nothing
// Errors: none
void debugPrint(const std::string& msg) {
    if (DEBUG) {
        std::cout << "[DEBUG] " << msg << std::endl;
    }
}
