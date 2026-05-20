#include "algorithm.h"
#include "helpers.h"
#include "memoization.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Keep debug message construction out of the recursive hot path when DEBUG is off.
#define debugPrint(message_expr) \
    do {                         \
        if (DEBUG) {             \
            ::debugPrint(message_expr); \
        }                        \
    } while (false)

namespace {

enum class HeuristicStateKind : std::uint8_t {
    Pair = 0,
    S = 1,
    I = 2,
};

struct CachedFreeEdge {
    bool has_edge = false;
    std::pair<int, int> edge{-1, -1};
};

struct HeuristicBundle {
    int conflict_lb = 0;
    int state_bounds_lb = 0;
    int small_exact_lb = 0;
    int local_nopath_lb = 0;
    int branchy_core_lb = 0;
    int combined_lb = 0;
};

struct EmptySCoreSummary {
    int lower_bound = 0;
    int exact_total = 0;
    bool fully_exact = false;
    int common_splits = 0;
    int free_reductions = 0;
    std::string motif = "unknown";
    std::string signature;
};

struct CoreDecompositionCandidate {
    VectorRangeTreeMap A1;
    VectorRangeTreeMap A2;
    VectorRangeTreeMap B1;
    VectorRangeTreeMap B2;
    int score = 0;
    int peak = 0;
    int min_side_edges = 0;
};

struct FreeEdgeCoreReduction {
    bool valid = false;
    bool solved = false;
    VectorRangeTreeMap left_start;
    VectorRangeTreeMap left_target;
    VectorRangeTreeMap right_start;
    VectorRangeTreeMap right_target;
};

struct TreeStructureProfile {
    int branching_nodes = 0;
    int max_branch_subtree_edges = 0;
    int chain_arm_count = 0;
    int broom_arm_count = 0;
    std::string degree_histogram = "0:0|1:0|2:0";
    std::string motif = "unknown";
};

struct MotifPeelLayout {
    VectorRangeTreeMap peeled_start;
    VectorRangeTreeMap peeled_target;
    VectorRangeTreeMap residual_start;
    VectorRangeTreeMap residual_target;
    int base_cost = 0;
    int peeled_edges = 0;
    int residual_edges = 0;
    bool solves_entire_state = false;
    std::string peel_motif = "unknown";
    std::string kind;
    std::string articulation_edge;
};

struct ArticulationArmReduction {
    bool valid = false;
    bool solved = false;
    VectorRangeTreeMap residual_start;
    VectorRangeTreeMap residual_target;
    int exact_side_cost = 0;
    int forced_cost = 0;
    int peeled_edges = 0;
    int residual_edges = 0;
    std::string residual_motif = "unknown";
    std::string articulation_edge;
    std::string peel_motif = "unknown";
};

struct EmptySMotifSummary {
    int reduced_core_internal_edges = 0;
    int chain_arm_count = 0;
    int broom_arm_count = 0;
    std::string motif = "unknown";
    std::string signature;
    TreeStructureProfile start_profile;
    TreeStructureProfile target_profile;
};

struct BranchyCoreSummary {
    bool reduction_triggered = false;
    int forced_cost = 0;
    int residual_edges = 0;
    int residual_conflicts = 0;
    int exact_bound = -1;
    bool exact_bound_used = false;
    std::string residual_motif = "unknown";
    std::string residual_signature;
};

struct PrefixPairRequestKey {
    Key128 tree_key{};
    std::uint64_t diagonal_key = 0;
    bool operator==(const PrefixPairRequestKey& o) const {
        return tree_key == o.tree_key && diagonal_key == o.diagonal_key;
    }
};

struct PrefixPairRequestKeyHash {
    std::size_t operator()(const PrefixPairRequestKey& k) const {
        std::size_t h1 = Key128Hash{}(k.tree_key);
        std::size_t h2 = std::hash<std::uint64_t>{}(k.diagonal_key);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

struct PrefixPairRequestMemoEntry {
    int count = 0;
    std::array<std::pair<std::pair<int, int>, std::pair<int, int>>, 2> pairs{};
    int vcount = 0;
    int bucket_size_l = 0;
    int bucket_size_r = 0;
};

thread_local std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> g_prefix_pair_scratch;

struct FreeEdgePartitionStructureKey {
    Key128 pair_key{};
    std::uint32_t parent_range_l = 0;
    std::uint32_t parent_range_r = 0;
    std::uint32_t child_range_l = 0;
    std::uint32_t child_range_r = 0;
    bool operator==(const FreeEdgePartitionStructureKey& o) const {
        return pair_key == o.pair_key &&
               parent_range_l == o.parent_range_l && parent_range_r == o.parent_range_r &&
               child_range_l == o.child_range_l && child_range_r == o.child_range_r;
    }
};

struct FreeEdgePartitionStructureKeyHash {
    std::size_t operator()(const FreeEdgePartitionStructureKey& k) const {
        std::size_t h = Key128Hash{}(k.pair_key);
        h ^= std::hash<std::uint32_t>{}(k.parent_range_l) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<std::uint32_t>{}(k.parent_range_r) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<std::uint32_t>{}(k.child_range_l) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<std::uint32_t>{}(k.child_range_r) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct FreeEdgePartitionStructureCacheEntry {
    bool computed = false;
    bool precheck_ok = true;
    bool partition_nodes_match = true;
    bool had_exception = false;
    int n1 = 0;
    int n2 = 0;
    VectorRangeTreeMap T_bar1;
    VectorRangeTreeMap T_bar2;
    VectorRangeTreeMap T_end1;
    VectorRangeTreeMap T_end2;
    const VectorRangeTreeMap* T_end1_ptr = nullptr;
    const VectorRangeTreeMap* T_end2_ptr = nullptr;
};

struct TargetPartitionRangeCacheKey {
    Key128 tree_key{};
    std::uint32_t child_range_l = 0;
    std::uint32_t child_range_r = 0;

    bool operator==(const TargetPartitionRangeCacheKey& o) const {
        return tree_key == o.tree_key &&
               child_range_l == o.child_range_l &&
               child_range_r == o.child_range_r;
    }
};

struct TargetPartitionRangeCacheKeyHash {
    std::size_t operator()(const TargetPartitionRangeCacheKey& k) const {
        std::size_t h = Key128Hash{}(k.tree_key);
        h ^= std::hash<std::uint32_t>{}(k.child_range_l) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<std::uint32_t>{}(k.child_range_r) + 0x517cc1b727220a95ULL + (h << 5) + (h >> 3);
        return h;
    }
};

struct TargetPartitionRangeCacheEntry {
    bool computed = false;
    VectorRangeTreeMap T_end1;
    VectorRangeTreeMap T_end2;
};

struct FreeEdgePartitionSideCacheKey {
    FreeEdgePartitionStructureKey structure_key{};
    Key128 bounds_key{};
    bool operator==(const FreeEdgePartitionSideCacheKey& o) const {
        return structure_key == o.structure_key && bounds_key == o.bounds_key;
    }
};

struct FreeEdgePartitionSideCacheKeyHash {
    std::size_t operator()(const FreeEdgePartitionSideCacheKey& k) const {
        std::size_t h = FreeEdgePartitionStructureKeyHash{}(k.structure_key);
        h ^= Key128Hash{}(k.bounds_key) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct FreeEdgePartitionSideCacheEntry {
    bool computed = false;
    bool had_exception = false;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> S1;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> S2;
    int conf1 = 0;
    int conf2 = 0;
    Key128 s1_pair_key{};
    Key128 s2_pair_key{};
    Key128 s1_bounds_key{};
    Key128 s2_bounds_key{};
};

struct EmptySChildSignature {
    Key128 child_pair_key{};
    Key128 child_bounds_key{};
    int next_conflicts_bucket = INT_MAX;
    bool operator==(const EmptySChildSignature& o) const {
        return child_pair_key == o.child_pair_key &&
               child_bounds_key == o.child_bounds_key &&
               next_conflicts_bucket == o.next_conflicts_bucket;
    }
};

struct EmptySChildSignatureHash {
    std::size_t operator()(const EmptySChildSignature& s) const {
        std::size_t h = Key128Hash{}(s.child_pair_key);
        h ^= Key128Hash{}(s.child_bounds_key) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(s.next_conflicts_bucket) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct EmptySChildCacheItem {
    std::pair<int, int> edge{};
    int next_conflicts = INT_MAX;
    int next_conflicts_bucket = INT_MAX;
    int free_hint = 0;
    int child_lb = INT_MAX;
    int child_k = -1;
    int pair_incumbent = INT_MAX;
    Key128 child_pair_key{};
    Key128 child_bounds_key{};
    std::uint64_t partition_signature = 0;
};

struct EmptySChildCacheEntry {
    std::vector<EmptySChildCacheItem> items;
};

struct EmptySChildCacheKey {
    Key128 pair_key{};

    bool operator==(const EmptySChildCacheKey& o) const {
        return pair_key == o.pair_key;
    }
};

struct EmptySChildCacheKeyHash {
    std::size_t operator()(const EmptySChildCacheKey& k) const {
        std::size_t h = Key128Hash{}(k.pair_key);
        return h;
    }
};

struct FreeEdgePartitionSplitSignatureCacheKey {
    Key128 s1_bounds_key{};
    Key128 s2_bounds_key{};

    bool operator==(const FreeEdgePartitionSplitSignatureCacheKey& o) const {
        return s1_bounds_key == o.s1_bounds_key &&
               s2_bounds_key == o.s2_bounds_key;
    }
};

struct FreeEdgePartitionSplitSignatureCacheKeyHash {
    std::size_t operator()(const FreeEdgePartitionSplitSignatureCacheKey& k) const {
        std::size_t h = Key128Hash{}(k.s1_bounds_key);
        h ^= Key128Hash{}(k.s2_bounds_key) + 0x517cc1b727220a95ULL + (h << 5) + (h >> 3);
        return h;
    }
};

bool key128Less(const Key128& a, const Key128& b) {
    if (a.hi != b.hi) return a.hi < b.hi;
    return a.lo < b.lo;
}

bool refreshOriginalPreorderInPlace(VectorRangeTreeMap& T) {
    if (T.root < 0 || T.original_nodes.empty() || !T.isOriginal(T.root)) {
        return false;
    }
    T.original_preorder.clear();
    T.original_preorder.reserve(T.original_nodes.size());

    auto dfs = [&](auto&& self, int node) -> void {
        if (node < 0 || !T.isOriginal(node)) return;
        T.original_preorder.push_back(node);
        const int left = T.getLeftChild(node);
        const int right = T.getRightChild(node);
        if (left >= 0 && T.isOriginal(left)) self(self, left);
        if (right >= 0 && T.isOriginal(right)) self(self, right);
    };
    dfs(dfs, T.root);
    return T.original_preorder.size() == T.original_nodes.size();
}

FreeEdgePartitionSplitSignatureCacheKey makeSplitSignatureCacheKey(const Key128& a,
                                                                   const Key128& b) {
    if (key128Less(b, a)) {
        return FreeEdgePartitionSplitSignatureCacheKey{b, a};
    }
    return FreeEdgePartitionSplitSignatureCacheKey{a, b};
}

struct FreeEdgePartitionSplitDecisionCacheKey {
    FreeEdgePartitionSideCacheKey side_cache_key{};
    int remaining_budget = -1;

    bool operator==(const FreeEdgePartitionSplitDecisionCacheKey& o) const {
        return side_cache_key == o.side_cache_key &&
               remaining_budget == o.remaining_budget;
    }
};

struct FreeEdgePartitionSplitDecisionCacheKeyHash {
        std::size_t operator()(const FreeEdgePartitionSplitDecisionCacheKey& k) const {
        std::size_t h = FreeEdgePartitionSideCacheKeyHash{}(k.side_cache_key);
        h ^= std::hash<int>{}(k.remaining_budget) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct FreeEdgePartitionSplitDecisionCacheEntry {
    bool computed = false;
    bool solvable = false;
    bool had_exception = false;
    int witness_k1 = -1;
};

struct FreeEdgePartitionSplitBudgetCacheEntry {
    bool computed = false;
    int solved_from = -1;
    int failed_up_to = -1;
    int witness_k1 = -1;
};

struct FreeEdgePartitionSideBudgetRangeCacheKey {
    Key128 side_pair_key{};
    Key128 side_bounds_key{};
    std::uint32_t range_lo = 0;
    std::uint32_t range_hi = 0;

    bool operator==(const FreeEdgePartitionSideBudgetRangeCacheKey& o) const {
        return side_pair_key == o.side_pair_key &&
               side_bounds_key == o.side_bounds_key &&
               range_lo == o.range_lo &&
               range_hi == o.range_hi;
    }
};

struct FreeEdgePartitionSideBudgetRangeCacheKeyHash {
    std::size_t operator()(const FreeEdgePartitionSideBudgetRangeCacheKey& k) const {
        std::size_t h = Key128Hash{}(k.side_pair_key);
        h ^= Key128Hash{}(k.side_bounds_key) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<std::uint32_t>{}(k.range_lo) + 0x517cc1b727220a95ULL + (h << 5) + (h >> 3);
        h ^= std::hash<std::uint32_t>{}(k.range_hi) + 0x517cc1b727220a95ULL + (h << 4) + (h >> 4);
        return h;
    }
};

struct FreeEdgePartitionSideBudgetRangeCacheEntry {
    bool computed = false;
    int min_feasible_budget = std::numeric_limits<int>::max();
};

struct FreeEdgePartitionSideBudgetCacheKey {
    Key128 side_pair_key{};
    Key128 side_bounds_key{};

    bool operator==(const FreeEdgePartitionSideBudgetCacheKey& o) const {
        return side_pair_key == o.side_pair_key &&
               side_bounds_key == o.side_bounds_key;
    }
};

struct FreeEdgePartitionSideBudgetCacheKeyHash {
    std::size_t operator()(const FreeEdgePartitionSideBudgetCacheKey& k) const {
        std::size_t h = Key128Hash{}(k.side_pair_key);
        h ^= Key128Hash{}(k.side_bounds_key) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct FreeEdgePartitionSideBudgetCacheEntry {
    bool computed = false;
    int known_not_feasible_below = -1;
    int known_min_feasible_budget = std::numeric_limits<int>::max();
};

struct CommonEdgeDecomposeCacheKey {
    Key128 pair_key{};
    int budget = 0;
    bool operator==(const CommonEdgeDecomposeCacheKey& o) const {
        return pair_key == o.pair_key && budget == o.budget;
    }
};

struct CommonEdgeDecomposeCacheKeyHash {
    std::size_t operator()(const CommonEdgeDecomposeCacheKey& k) const {
        std::size_t h = Key128Hash{}(k.pair_key);
        h ^= std::hash<int>{}(k.budget) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct CommonEdgeDecomposeCacheValue {
    bool solvable = false;
    bool handled = false;
};

std::unordered_map<Key128, int, Key128Hash> g_conflict_cache;
std::unordered_map<Key128, std::unordered_set<RP, RPH>, Key128Hash> g_target_edge_set_cache;
std::unordered_map<Key128, CachedFreeEdge, Key128Hash> g_free_edge_cache;
std::unordered_map<Key128, HeuristicBundle, Key128Hash> g_heuristic_bundle_cache;
std::unordered_map<Key128, int, Key128Hash> g_small_exact_dist_cache;
std::unordered_map<Key128, int, Key128Hash> g_local_nopath_lb_cache;
std::unordered_map<Key128, int, Key128Hash> g_tds_inflight_best_budget;
std::unordered_map<Key128, int, Key128Hash> g_tdi_inflight_best_budget;
std::unordered_map<Key128, int, Key128Hash> g_pair_incumbent_upper;
std::size_t g_active_problem_size = 0;
std::unordered_map<Key128, EmptySCoreSummary, Key128Hash> g_empty_s_core_cache;
std::unordered_map<Key128, int, Key128Hash> g_empty_s_core_exact_cache;
std::unordered_map<Key128, EmptySMotifSummary, Key128Hash> g_empty_s_motif_summary_cache;
std::unordered_map<Key128, int, Key128Hash> g_empty_s_motif_exact_cache;
std::unordered_map<Key128, BranchyCoreSummary, Key128Hash> g_branchy_core_summary_cache;
std::unordered_map<Key128, int, Key128Hash> g_branchy_core_exact_cache;
std::unordered_map<Key128, bool, Key128Hash> g_pair_common_edge_cache;
std::unordered_map<Key128, bool, Key128Hash> g_empty_s_free_hint_gate_cache;
std::unordered_map<PrefixPairRequestKey, PrefixPairRequestMemoEntry, PrefixPairRequestKeyHash> g_prefix_pair_request_memo;
std::unordered_map<EmptySChildCacheKey, EmptySChildCacheEntry, EmptySChildCacheKeyHash> g_empty_s_child_cache;
std::unordered_map<FreeEdgePartitionStructureKey, FreeEdgePartitionStructureCacheEntry, FreeEdgePartitionStructureKeyHash>
    g_free_edge_partition_structure_cache;
std::unordered_map<TargetPartitionRangeCacheKey, TargetPartitionRangeCacheEntry, TargetPartitionRangeCacheKeyHash>
    g_target_partition_range_cache;
std::unordered_map<FreeEdgePartitionSideCacheKey, FreeEdgePartitionSideCacheEntry, FreeEdgePartitionSideCacheKeyHash>
    g_free_edge_partition_side_cache;
std::unordered_map<FreeEdgePartitionSplitDecisionCacheKey,
                   FreeEdgePartitionSplitDecisionCacheEntry,
                   FreeEdgePartitionSplitDecisionCacheKeyHash>
    g_free_edge_partition_split_cache;
std::unordered_map<FreeEdgePartitionSplitSignatureCacheKey,
                   FreeEdgePartitionSplitBudgetCacheEntry,
                   FreeEdgePartitionSplitSignatureCacheKeyHash>
    g_free_edge_partition_split_signature_cache;
std::unordered_map<FreeEdgePartitionSideCacheKey, FreeEdgePartitionSplitBudgetCacheEntry, FreeEdgePartitionSideCacheKeyHash>
    g_free_edge_partition_split_budget_cache;
std::unordered_map<FreeEdgePartitionSideBudgetCacheKey,
                   FreeEdgePartitionSideBudgetCacheEntry,
                   FreeEdgePartitionSideBudgetCacheKeyHash>
    g_free_edge_partition_side_budget_cache;
std::unordered_map<FreeEdgePartitionSideBudgetRangeCacheKey,
                   FreeEdgePartitionSideBudgetRangeCacheEntry,
                   FreeEdgePartitionSideBudgetRangeCacheKeyHash>
    g_free_edge_partition_side_budget_range_cache;
std::unordered_map<CommonEdgeDecomposeCacheKey, CommonEdgeDecomposeCacheValue, CommonEdgeDecomposeCacheKeyHash>
    g_common_edge_decompose_cache;
std::unordered_map<CommonEdgeDecomposeCacheKey, bool, CommonEdgeDecomposeCacheKeyHash>
    g_common_edge_double_rotate_cache;
std::unordered_set<Key128, Key128Hash> g_plateau_tree_dumped_states;
int g_plateau_tree_dump_count = 0;
thread_local int g_empty_s_disable_free_hint_depth = 0;

bool branchyCoreLbEnabled();
int branchyCoreExactM();
int smallExactBoundM();
int localNoPathDepth();
BranchyCoreSummary analyzeBranchyCoreSummary(const VectorRangeTreeMap& start,
                                             const VectorRangeTreeMap& target,
                                             int exact_m);

bool articulationArmReduceEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_ARTIC_ARM_REDUCE");
    cached = (env && std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

int articulationArmMaxEdges() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_ARTIC_ARM_MAX_EDGES");
    if (!env || std::string(env).empty()) {
        cached = 5;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 0) value = 0;
    if (value > 10) value = 10;
    cached = value;
    return cached;
}

const std::string& articulationArmLogFilePath() {
    static std::string cached = []() {
        const char* env = std::getenv("FLIPDIST_ARTIC_ARM_LOG_FILE");
        return env ? std::string(env) : std::string();
    }();
    return cached;
}

bool articulationArmLogEnabled() {
    return !articulationArmLogFilePath().empty();
}

bool tdiPrefixPruneEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_TDI_PREFIX_PRUNE");
    cached = (env && std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

int countCommonEdges(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end,
                     int stop_after = 1) {
    if (T_init.original_nodes.empty() || T_end.original_nodes.empty()) {
        return 0;
    }
    int count = 0;
    try {
        const auto target_edges = buildTargetSet(T_end);
        const int threshold = std::max(1, stop_after);
        for (int v : T_init.original_nodes) {
            if (!T_init.isOriginal(v)) {
                continue;
            }
            auto pr = T_init.getRange(v);
            int left = T_init.getLeftChild(v);
            if (left >= 0 && T_init.isOriginal(left)) {
                auto cr = T_init.getRange(left);
                RP edge{pr.first, pr.second, cr.first, cr.second};
                if (target_edges.count(edge)) {
                    if (++count >= threshold) {
                        return count;
                    }
                }
            }
            int right = T_init.getRightChild(v);
            if (right >= 0 && T_init.isOriginal(right)) {
                auto cr = T_init.getRange(right);
                RP edge{pr.first, pr.second, cr.first, cr.second};
                if (target_edges.count(edge)) {
                    if (++count >= threshold) {
                        return count;
                    }
                }
            }
        }
    } catch (const SearchAbortException&) {
        throw;
    } catch (...) {
        return 0;
    }
    return count;
}

int tdsCommonEdgeMinCount() {
    static int cached = -1;
    if (cached >= 0) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_TDS_COMMON_EDGE_MIN");
    if (!env || std::string(env).empty()) {
        cached = 3;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 1) {
        value = 1;
    }
    cached = value;
    return cached;
}

bool pairHasCommonEdgesCached(const VectorRangeTreeMap& T_a,
                              const VectorRangeTreeMap& T_b,
                              int threshold) {
    if (threshold <= 0) return true;
    const Key128 pair_key = makeKeyPair(T_a, T_b);
    auto it = g_pair_common_edge_cache.find(pair_key);
    if (it != g_pair_common_edge_cache.end()) {
        return it->second;
    }
    const bool has_common = (countCommonEdges(T_a, T_b, threshold) >= threshold);
    g_pair_common_edge_cache[pair_key] = has_common;
    return has_common;
}

int tdiPrefixPruneMinSSize() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_TDI_PREFIX_MIN_S");
    if (!env || std::string(env).empty()) {
        cached = 2;
    } else {
        cached = std::max(1, std::atoi(env));
    }
    return cached;
}

int cachedConflictCount(const VectorRangeTreeMap& start, const VectorRangeTreeMap& target);
int cachedConflictCountWithKey(const VectorRangeTreeMap& start,
                               const VectorRangeTreeMap& target,
                               const Key128& key);
int pairIncumbentUpper(const Key128& pair_key);

bool tdsPartitionCacheEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* disable_env = std::getenv("FLIPDIST_DEBUG_DISABLE_PARTITION_CACHE");
    if (disable_env && std::string(disable_env) == "1") {
        cached = 0;
    } else if (const char* enable_env = std::getenv("FLIPDIST_DEBUG_ENABLE_PARTITION_CACHE");
               enable_env && std::string(enable_env) == "1") {
        cached = 1;
    } else {
        cached = 0;
    }
    return cached == 1;
}

bool tdsPartitionSplitCacheEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* disable_env = std::getenv("FLIPDIST_DEBUG_DISABLE_PARTITION_SPLIT_CACHE");
    if (disable_env && std::string(disable_env) == "1") {
        cached = 0;
    } else if (const char* enable_env = std::getenv("FLIPDIST_DEBUG_ENABLE_PARTITION_SPLIT_CACHE");
               enable_env && std::string(enable_env) == "1") {
        cached = 1;
    } else {
        cached = 0;
    }
    return cached == 1;
}

bool tdsPartitionMinBudgetSearchEnabledForSize(const std::size_t n) {
    static int cached = -1;
    if (cached >= 0) {
        if (cached == 0) return false;
        if (cached == 1) return true;
        return n >= 24;
    }
    const char* env = std::getenv("FLIPDIST_TDS_PARTITION_MIN_BUDGET_SEARCH");
    if (env && std::string(env) == "0") {
        cached = 0;
    } else if (env && std::string(env) == "1") {
        cached = 1;
    } else {
        cached = 2;
    }
    if (cached == 0) return false;
    if (cached == 1) return true;
    return n >= 24;
}

bool tdsPartitionLowerFirstEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* env = std::getenv("FLIPDIST_TDS_PARTITION_LOWER_FIRST");
    if (env && std::string(env) == "0") {
        cached = 0;
    } else {
        cached = 1;
    }
    return cached == 1;
}

int tdsPartitionMinSideMode() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) return cached;
    const char* env = std::getenv("FLIPDIST_TDS_PARTITION_MIN_SIDE");
    const std::string mode = env ? std::string(env) : std::string();
    if (mode == "1" || mode == "side1") {
        cached = 1;
    } else if (mode == "2" || mode == "side2") {
        cached = 2;
    } else if (mode == "lb" || mode == "3") {
        cached = 3;
    } else if (mode == "conflict" || mode == "4") {
        cached = 4;
    } else if (mode == "size" || mode == "5") {
        cached = 5;
    } else if (mode == "pairs" || mode == "6") {
        cached = 6;
    } else {
        cached = 1;
    }
    return cached;
}

bool emptySDedupEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* env = std::getenv("FLIPDIST_DEBUG_DISABLE_EMPTY_S_DEDUPE");
    if (env && std::string(env) == "1") {
        cached = 0;
    } else {
        cached = 1;
    }
    return cached == 1;
}

int emptySOrderMode() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) return cached;
    const char* env = std::getenv("FLIPDIST_EMPTY_S_ORDER");
    const std::string mode = env ? std::string(env) : std::string();
    if (mode.empty() || mode == "lb") {
        cached = 0;
    } else if (mode == "free") {
        cached = 1;
    } else if (mode == "conflict") {
        cached = 2;
    } else {
        cached = 0;
    }
    return cached;
}

bool emptySTightNoCacheEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* env = std::getenv("FLIPDIST_EMPTY_S_TIGHT_NO_CACHE");
    if (env && std::string(env) == "0") {
        cached = 0;
    } else {
        cached = 1;
    }
    return cached == 1;
}

bool emptySRootShapeFreeHintGateEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* env = std::getenv("FLIPDIST_EMPTY_S_ROOT_SHAPE_FREE_HINT_GATE");
    if (env && std::string(env) == "0") {
        cached = 0;
    } else {
        cached = 1;
    }
    return cached == 1;
}

struct EmptySShapeStats {
    int height = 0;
    int branching = 0;
    double avg_depth = 0.0;
};

EmptySShapeStats emptySShapeStats(const VectorRangeTreeMap& T) {
    EmptySShapeStats stats;
    long long depth_sum = 0;
    int node_count = 0;
    auto dfs = [&](auto&& self, int node, int depth) -> void {
        if (node < 0 || !T.isOriginal(node)) {
            return;
        }
        ++node_count;
        depth_sum += depth;
        stats.height = std::max(stats.height, depth);
        const int left = T.getLeftChild(node);
        const int right = T.getRightChild(node);
        if (left >= 0 && T.isOriginal(left) && right >= 0 && T.isOriginal(right)) {
            ++stats.branching;
        }
        self(self, left, depth + 1);
        self(self, right, depth + 1);
    };
    dfs(dfs, T.root, 0);
    if (node_count > 0) {
        stats.avg_depth = static_cast<double>(depth_sum) / static_cast<double>(node_count);
    }
    return stats;
}

bool shouldDisableEmptySFreeHintsForRootShape(const VectorRangeTreeMap& start,
                                              const VectorRangeTreeMap& target,
                                              const Key128& pair_key) {
    if (!emptySRootShapeFreeHintGateEnabled()) {
        return false;
    }
    if (start.original_nodes.size() < 25 || target.original_nodes.size() < 25) {
        return false;
    }
    auto cache_it = g_empty_s_free_hint_gate_cache.find(pair_key);
    if (cache_it != g_empty_s_free_hint_gate_cache.end()) {
        return cache_it->second;
    }
    const EmptySShapeStats source = emptySShapeStats(start);
    const EmptySShapeStats dest = emptySShapeStats(target);
    const bool should_disable =
        source.branching >= dest.branching + 2 &&
        source.avg_depth + 0.1 < dest.avg_depth &&
        dest.height <= source.height;
    g_empty_s_free_hint_gate_cache.emplace(pair_key, should_disable);
    return should_disable;
}

struct EmptySFreeHintDisableGuard {
    bool active = false;
    explicit EmptySFreeHintDisableGuard(bool enabled) : active(enabled) {
        if (active) {
            ++g_empty_s_disable_free_hint_depth;
        }
    }
    ~EmptySFreeHintDisableGuard() {
        if (active) {
            --g_empty_s_disable_free_hint_depth;
        }
    }
};

bool tdsInflightPruningEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* disable_env = std::getenv("FLIPDIST_DEBUG_DISABLE_TDS_INFLIGHT");
    if (disable_env && std::string(disable_env) == "1") {
        cached = 0;
        return cached == 1;
    }
    const char* enable_env = std::getenv("FLIPDIST_ENABLE_TDS_INFLIGHT");
    if (enable_env && std::string(enable_env) == "1") {
        cached = 1;
        return cached == 1;
    }
    if (enable_env && std::string(enable_env) == "0") {
        cached = 0;
        return cached == 1;
    }
    if (!enable_env) {
        cached = 1;
    } else {
        cached = 1;
    }
    return cached == 1;
}

bool tdiInflightPruningEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* disable_env = std::getenv("FLIPDIST_DEBUG_DISABLE_TDI_INFLIGHT");
    if (disable_env && std::string(disable_env) == "1") {
        cached = 0;
        return cached == 1;
    }
    const char* enable_env = std::getenv("FLIPDIST_ENABLE_TDI_INFLIGHT");
    if (enable_env && std::string(enable_env) == "1") {
        cached = 1;
    } else {
        cached = 0;
    }
    return cached == 1;
}

bool tdsCommonEdgeRotatePairEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* disable_env = std::getenv("FLIPDIST_DEBUG_DISABLE_TDS_COMMON_EDGE_PAIR");
    if (disable_env && std::string(disable_env) == "1") {
        cached = 0;
        return cached == 1;
    }
    const char* enable_env = std::getenv("FLIPDIST_ENABLE_TDS_COMMON_EDGE_PAIR");
    if (enable_env && std::string(enable_env) == "1") {
        cached = 1;
    } else {
        cached = 0;
    }
    return cached == 1;
}

static std::vector<RP> collectCommonEdgesByTarget(
    const VectorRangeTreeMap& tree,
    const std::unordered_set<RP, RPH>& target_edges) {
    std::vector<RP> edges;
    try {
        if (tree.original_nodes.empty() || target_edges.empty()) {
            return edges;
        }
        edges.reserve(tree.original_nodes.size());
        for (int v : tree.original_nodes) {
            if (!tree.isOriginal(v)) {
                continue;
            }
            auto pr = tree.getRange(v);
            int left = tree.getLeftChild(v);
            if (left >= 0 && tree.isOriginal(left)) {
                auto cr = tree.getRange(left);
                RP edge{pr.first, pr.second, cr.first, cr.second};
                if (target_edges.count(edge)) {
                    edges.push_back(edge);
                }
            }
            int right = tree.getRightChild(v);
            if (right >= 0 && tree.isOriginal(right)) {
                auto cr = tree.getRange(right);
                RP edge{pr.first, pr.second, cr.first, cr.second};
                if (target_edges.count(edge)) {
                    edges.push_back(edge);
                }
            }
        }
    } catch (const SearchAbortException&) {
        throw;
    } catch (...) {
        edges.clear();
    }
    return edges;
}

bool tryCommonEdgeDecomposeAfterOneRotationEach(const VectorRangeTreeMap& T_start,
                                               const VectorRangeTreeMap& T_end,
                                               int k) {
    if (k <= 0 || T_start.original_nodes.empty() || T_end.original_nodes.empty()) {
        return false;
    }

    const int remaining_budget = k - 1;
    const int remaining_budget_after_both = k - 2;

    const auto start_edges = getInternalEdges(T_start);
    const auto end_edges = getInternalEdges(T_end);
    const int max_edges_start = static_cast<int>(start_edges.size());
    const int max_edges_end = static_cast<int>(end_edges.size());
    if (max_edges_start <= 0 || max_edges_end <= 0 || remaining_budget < 0) {
        return false;
    }

    const auto start_target_edges = buildTargetSet(T_end);
    const auto end_target_edges = buildTargetSet(T_start);

    struct RotatedWithCommonEdges {
        VectorRangeTreeMap tree;
        std::vector<RP> common_edges;
        Key128 pair_key{};
        int pair_lb = INT_MAX;
    };

    std::vector<RotatedWithCommonEdges> start_neighbors;
    std::vector<RotatedWithCommonEdges> end_neighbors;
    start_neighbors.reserve(max_edges_start);
    end_neighbors.reserve(max_edges_end);
    std::unordered_set<Key128, Key128Hash> seen_start_states;
    std::unordered_set<Key128, Key128Hash> seen_end_states;
    seen_start_states.reserve(max_edges_start * 2 + 1);
    seen_end_states.reserve(max_edges_end * 2 + 1);

    auto rotateEdgeOnce = [&](VectorRangeTreeMap& probe, int parent, int child, int& rotate_dir) -> bool {
        if (probe.getLeftChild(parent) == child) {
            probe.rotateRight(parent);
            rotate_dir = 1;
            return true;
        }
        if (probe.getRightChild(parent) == child) {
            probe.rotateLeft(parent);
            rotate_dir = -1;
            return true;
        }
        rotate_dir = 0;
        return false;
    };

    auto undoRotateState = [&](VectorRangeTreeMap& probe, int child, int rotate_dir) {
        if (rotate_dir == 1) {
            probe.rotateLeft(child);
        } else if (rotate_dir == -1) {
            probe.rotateRight(child);
        }
    };

    VectorRangeTreeMap start_probe = safeCopyTree(T_start);
    for (const auto& edge : start_edges) {
        int parent = edge.first;
        int child = edge.second;
        int rotate_dir = 0;
        if (!rotateEdgeOnce(start_probe, parent, child, rotate_dir)) {
            continue;
        }
        auto rollbackStart = [&]() {
            if (rotate_dir != 0) {
                undoRotateState(start_probe, child, rotate_dir);
            }
        };
        try {
            Key128 pair_key = makeKeyPair(start_probe, T_end);
            if (pairIncumbentUpper(pair_key) <= remaining_budget) {
                rollbackStart();
                return true;
            }
            int pair_lb = std::max(cachedConflictCount(start_probe, T_end),
                                   requiredBudgetFromBounds(g_kbounds, pair_key));
            if (pair_lb > remaining_budget) {
                rollbackStart();
                continue;
            }
            if (!seen_start_states.insert(pair_key).second) {
                rollbackStart();
                continue;
            }
            auto common_edges = collectCommonEdgesByTarget(start_probe, start_target_edges);
            if (!common_edges.empty()) {
                start_neighbors.push_back(
                    RotatedWithCommonEdges{std::move(start_probe), std::move(common_edges), pair_key, pair_lb});
                start_probe = safeCopyTree(T_start);
            } else {
                rollbackStart();
            }
        } catch (const SearchAbortException&) {
            throw;
        } catch (...) {
            rollbackStart();
            continue;
        }
    }

    VectorRangeTreeMap end_probe = safeCopyTree(T_end);
    for (const auto& edge : end_edges) {
        int parent = edge.first;
        int child = edge.second;
        int rotate_dir = 0;
        if (!rotateEdgeOnce(end_probe, parent, child, rotate_dir)) {
            continue;
        }
        auto rollbackEnd = [&]() {
            if (rotate_dir != 0) {
                undoRotateState(end_probe, child, rotate_dir);
            }
        };
        try {
            Key128 pair_key = makeKeyPair(T_start, end_probe);
            if (pairIncumbentUpper(pair_key) <= remaining_budget) {
                rollbackEnd();
                return true;
            }
            int pair_lb = std::max(cachedConflictCount(T_start, end_probe),
                                   requiredBudgetFromBounds(g_kbounds, pair_key));
            if (pair_lb > remaining_budget) {
                rollbackEnd();
                continue;
            }
            if (!seen_end_states.insert(pair_key).second) {
                rollbackEnd();
                continue;
            }
            auto common_edges = collectCommonEdgesByTarget(end_probe, end_target_edges);
            if (!common_edges.empty()) {
                end_neighbors.push_back(
                    RotatedWithCommonEdges{std::move(end_probe), std::move(common_edges), pair_key, pair_lb});
                end_probe = safeCopyTree(T_end);
            } else {
                rollbackEnd();
            }
        } catch (const SearchAbortException&) {
            throw;
        } catch (...) {
            rollbackEnd();
            continue;
        }
    }

    if (remaining_budget_after_both >= 0 && !start_neighbors.empty() && !end_neighbors.empty()) {
        std::unordered_set<Key128, Key128Hash> seen_dual_pairs;
        seen_dual_pairs.reserve(
            std::max(1, static_cast<int>(start_neighbors.size() * end_neighbors.size())));
        std::unordered_map<RP, std::vector<int>, RPH> end_edge_to_indices;
        end_edge_to_indices.reserve(std::max(1, static_cast<int>(end_neighbors.size()) * 2));
        for (int end_idx = 0; end_idx < static_cast<int>(end_neighbors.size()); ++end_idx) {
            auto& signatures = end_neighbors[end_idx].common_edges;
            for (const auto& edge_key : signatures) {
                end_edge_to_indices[edge_key].push_back(end_idx);
            }
        }

        std::vector<int> seen_end_indices(end_neighbors.size(), -1);
        int mark = 0;
        for (int start_idx = 0; start_idx < static_cast<int>(start_neighbors.size()); ++start_idx) {
            ++mark;
            auto& start_signature = start_neighbors[start_idx].common_edges;
            for (const auto& edge_key : start_signature) {
                auto it = end_edge_to_indices.find(edge_key);
                if (it == end_edge_to_indices.end()) {
                    continue;
                }
                for (int end_idx : it->second) {
                    if (end_idx < 0 || end_idx >= static_cast<int>(end_neighbors.size())) {
                        continue;
                    }
                    if (seen_end_indices[end_idx] == mark) {
                        continue;
                    }
                    seen_end_indices[end_idx] = mark;
                    const auto& start_state = start_neighbors[start_idx];
                    const auto& end_state = end_neighbors[end_idx];
                    const auto dual_key = makeKeyPair(start_state.tree, end_state.tree);
                    if (!seen_dual_pairs.insert(dual_key).second) {
                        continue;
                    }
                    if (pairIncumbentUpper(dual_key) <= remaining_budget_after_both) {
                        return true;
                    }
                    int dual_lb = std::max(cachedConflictCount(start_state.tree, end_state.tree),
                                          requiredBudgetFromBounds(g_kbounds, dual_key));
                    if (dual_lb > remaining_budget_after_both) {
                        continue;
                    }
                    bool handled = false;
                    if (start_state.pair_lb > remaining_budget || end_state.pair_lb > remaining_budget) {
                        continue;
                    }
                    if (tryCommonEdgeDecomposition(start_state.tree, end_state.tree, remaining_budget_after_both, handled)) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

bool tryCommonEdgeDecomposeCached(const VectorRangeTreeMap& T_start,
                                 const VectorRangeTreeMap& T_end,
                                 int k,
                                 bool& handled) {
    const Key128 pair_key = makeKeyPair(T_start, T_end);
    const CommonEdgeDecomposeCacheKey cache_key{pair_key, k};
    auto it = g_common_edge_decompose_cache.find(cache_key);
    if (it != g_common_edge_decompose_cache.end()) {
        handled = it->second.handled;
        return it->second.solvable;
    }

    bool local_handled = false;
    const bool solvable = tryCommonEdgeDecomposition(T_start, T_end, k, local_handled);
    g_common_edge_decompose_cache.emplace(cache_key, CommonEdgeDecomposeCacheValue{solvable, local_handled});
    handled = local_handled;
    return solvable;
}

bool tryCommonEdgeDecomposeAfterOneRotationEachCached(const VectorRangeTreeMap& T_start,
                                                     const VectorRangeTreeMap& T_end,
                                                     int k) {
    const Key128 pair_key = makeKeyPair(T_start, T_end);
    const CommonEdgeDecomposeCacheKey cache_key{pair_key, k};
    auto it = g_common_edge_double_rotate_cache.find(cache_key);
    if (it != g_common_edge_double_rotate_cache.end()) {
        return it->second;
    }
    const bool solvable = tryCommonEdgeDecomposeAfterOneRotationEach(T_start, T_end, k);
    g_common_edge_double_rotate_cache.emplace(cache_key, solvable);
    return solvable;
}

std::string jsonEscape(const std::string& text) {
    std::string out;
    out.reserve(text.size() + 8);
    for (unsigned char c : text) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    std::ostringstream hex;
                    hex << "\\u" << std::hex << std::uppercase;
                    hex.width(4);
                    hex.fill('0');
                    hex << static_cast<unsigned int>(c);
                    out += hex.str();
                } else {
                    out.push_back(static_cast<char>(c));
                }
                break;
        }
    }
    return out;
}

std::pair<int, int> canonicalSPairEdge(std::pair<int, int> edge) {
    if (edge.second < edge.first) {
        std::swap(edge.first, edge.second);
    }
    return edge;
}

std::pair<std::pair<int, int>, std::pair<int, int>> canonicalSPair(
    std::pair<std::pair<int, int>, std::pair<int, int>> pair) {
    pair.first = canonicalSPairEdge(pair.first);
    pair.second = canonicalSPairEdge(pair.second);
    if (pair.second < pair.first) {
        std::swap(pair.first, pair.second);
    }
    return pair;
}

Key128 canonicalSPairKey(const std::pair<std::pair<int, int>, std::pair<int, int>>& pair) {
    auto canonical = canonicalSPair(pair);
    const std::uint64_t k1 =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(canonical.first.first)) << 32) ^
        static_cast<std::uint32_t>(canonical.first.second);
    const std::uint64_t k2 =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(canonical.second.first)) << 32) ^
        static_cast<std::uint32_t>(canonical.second.second);
    return Key128{k1, k2};
}

void appendUniquePartnerPairsFromDiagonal(
    const VectorRangeTreeMap& tree,
    const std::pair<int, int>& diagonal,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& out_pairs,
    std::unordered_set<Key128, Key128Hash>& seen_pairs,
    bool record_prefix_profile = true) {
    const Key128 tree_key = makeKeyPair(tree, tree);
    const int diag_a = std::min(diagonal.first, diagonal.second);
    const int diag_b = std::max(diagonal.first, diagonal.second);
    const std::uint64_t diagonal_key =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(diag_a)) << 32) ^
        static_cast<std::uint32_t>(diag_b);
    PrefixPairRequestKey request_key{tree_key, diagonal_key};
    PrefixPairgenRequestStats request_stats;
    int pair_count = 0;
    std::array<std::pair<std::pair<int, int>, std::pair<int, int>>, 2> pairs{};

    const auto memo_it = g_prefix_pair_request_memo.find(request_key);
    const bool memo_hit = (memo_it != g_prefix_pair_request_memo.end());
    if (!memo_hit) {
        g_prefix_pair_scratch.clear();
        auto helper_t0 = g_profile.enabled ? std::chrono::steady_clock::now()
                                           : std::chrono::steady_clock::time_point{};
        appendPartnerPairsFromSingleDiagonalProfiled(
            tree,
            diagonal,
            g_prefix_pair_scratch,
            &request_stats);
        if (record_prefix_profile && g_profile.enabled) {
            g_profile.time_tdi_prefix_pairgen_helper_ms +=
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - helper_t0).count();
        }
        pair_count = std::min<int>(2, static_cast<int>(g_prefix_pair_scratch.size()));
        for (int i = 0; i < pair_count; ++i) {
            pairs[i] = canonicalSPair(g_prefix_pair_scratch[i]);
        }
        PrefixPairRequestMemoEntry memo_entry;
        memo_entry.count = pair_count;
        memo_entry.pairs = pairs;
        memo_entry.vcount = request_stats.vcount;
        memo_entry.bucket_size_l = request_stats.bucket_size_l;
        memo_entry.bucket_size_r = request_stats.bucket_size_r;
        g_prefix_pair_request_memo.emplace(request_key, memo_entry);
        request_stats.emitted_pairs = pair_count;
    } else {
        const PrefixPairRequestMemoEntry& memo_entry = memo_it->second;
        pair_count = memo_entry.count;
        pairs = memo_entry.pairs;
        request_stats.emitted_pairs = memo_entry.count;
        request_stats.vcount = memo_entry.vcount;
        request_stats.bucket_size_l = memo_entry.bucket_size_l;
        request_stats.bucket_size_r = memo_entry.bucket_size_r;
    }

    if (record_prefix_profile && g_profile.enabled) {
        g_profile.tdi_prefix_request_total++;
        if (memo_hit) {
            g_profile.tdi_prefix_request_repeated++;
        } else {
            g_profile.tdi_prefix_request_unique++;
        }
        g_profile.tdi_prefix_request_vcount_samples.push_back(request_stats.vcount);
        g_profile.tdi_prefix_request_bucket_l_samples.push_back(request_stats.bucket_size_l);
        g_profile.tdi_prefix_request_bucket_r_samples.push_back(request_stats.bucket_size_r);
        if (request_stats.emitted_pairs <= 0) {
            g_profile.tdi_prefix_emitted_pairs_0++;
        } else if (request_stats.emitted_pairs == 1) {
            g_profile.tdi_prefix_emitted_pairs_1++;
        } else {
            g_profile.tdi_prefix_emitted_pairs_2++;
        }
        g_profile.tdi_prefix_pairgen_pairs_scanned += static_cast<long long>(pair_count);
    }
    auto insert_t0 = g_profile.enabled ? std::chrono::steady_clock::now()
                                       : std::chrono::steady_clock::time_point{};
    long long inserted = 0;
    for (int i = 0; i < pair_count; ++i) {
        auto& canonical = pairs[i];
        Key128 key = canonicalSPairKey(canonical);
        if (seen_pairs.insert(key).second) {
            out_pairs.push_back(canonical);
            inserted++;
        }
    }
    if (record_prefix_profile && g_profile.enabled) {
        g_profile.time_tdi_prefix_pairgen_insert_ms +=
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - insert_t0).count();
        g_profile.tdi_prefix_pairgen_pairs_inserted += inserted;
    }
}

void appendArticulationArmLogLine(const std::string& line) {
    const std::string& path = articulationArmLogFilePath();
    if (path.empty()) {
        return;
    }
    std::ofstream out(path, std::ios::app);
    if (!out) {
        return;
    }
    out << line << "\n";
}

Key128 mixHeuristicTag(Key128 key, std::uint64_t tag) {
    key.hi ^= (tag + 0x9e3779b97f4a7c15ULL + (key.hi << 6) + (key.hi >> 2));
    key.lo ^= ((tag << 1) + 0x517cc1b727220a95ULL + (key.lo << 5) + (key.lo >> 3));
    return key;
}

std::uint64_t heuristicConfigTag() {
    return (static_cast<std::uint64_t>(smallExactBoundM()) << 40) ^
           (static_cast<std::uint64_t>(localNoPathDepth()) << 32) ^
           (static_cast<std::uint64_t>(branchyCoreExactM()) << 1) ^
           static_cast<std::uint64_t>(branchyCoreLbEnabled() ? 1 : 0);
}

Key128 makePartitionFilteredSignature(Key128 pair_key, Key128 s_filtered_key) {
    Key128 key = mixHeuristicTag(pair_key, 0x7061727469746564ULL); // "partited"
    key.hi ^= (s_filtered_key.hi + 0x9e3779b97f4a7c15ULL + (key.hi << 6) + (key.hi >> 2));
    key.lo ^= (s_filtered_key.lo + 0x517cc1b727220a95ULL + (key.lo << 5) + (key.lo >> 3));
    return key;
}

int smallExactBoundM() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char *env = std::getenv("FLIPDIST_LB_SMALL_EXACT_M");
    if (!env || std::string(env).empty()) {
        cached = 0;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 0) value = 0;
    if (value > 30) value = 30;
    cached = value;
    return cached;
}

bool tdsPartitionPrecheckEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_TDS_PARTITION_PRECHECK");
    if (env && std::string(env) == "0") {
        cached = 0;
    } else {
        cached = 1;
    }
    return cached == 1;
}

int localNoPathDepth() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char *env = std::getenv("FLIPDIST_LB_LOCAL_DEPTH");
    if (!env || std::string(env).empty()) {
        cached = 0;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 0) value = 0;
    if (value > 2) value = 2;
    cached = value;
    return cached;
}

bool incumbentPruneEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_INCUMBENT_PRUNE");
    if (env && std::string(env) == "1") {
        cached = 1;
    } else {
        cached = 0;
    }
    return cached == 1;
}

bool combinedLowerBoundEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_DISABLE_COMBINED_LB");
    if (env && std::string(env) == "1") {
        cached = 0;
    } else {
        cached = 1;
    }
    return cached == 1;
}

bool defaultEmptySBoundConfig() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    cached = (smallExactBoundM() == 0 &&
              localNoPathDepth() == 0 &&
              !branchyCoreLbEnabled())
                 ? 1
                 : 0;
    return cached == 1;
}

int emptySLowerBoundFromKeys(const int conflict_lb,
                             const Key128& pair_key,
                             const Key128& bounds_key) {
    if (!combinedLowerBoundEnabled()) {
        return conflict_lb;
    }
    const int pair_lb = requiredBudgetFromBounds(g_kbounds, pair_key);
    const int state_lb = requiredBudgetFromBounds(g_tds_bounds, bounds_key);
    return std::max({conflict_lb, pair_lb, state_lb});
}

int minKProbeStep(const VectorRangeTreeMap& tree) {
    static int env_cached = INT_MIN;
    if (env_cached != INT_MIN) {
        return env_cached;
    }
    const char* env = std::getenv("FLIPDIST_MIN_K_PROBE_STEP");
    if (env && *env) {
        int value = std::atoi(env);
        if (value < 1) value = 1;
        if (value > 64) value = 64;
        env_cached = value;
        return env_cached;
    }
    const int n = static_cast<int>(tree.original_nodes.size());
    return std::max(4, (n + 1) / 4);
}

int minKLinearPrefixProbes(const VectorRangeTreeMap& tree,
                           const VectorRangeTreeMap& target) {
    static int env_cached = INT_MIN;
    if (env_cached != INT_MIN) {
        return env_cached;
    }
    const char* env = std::getenv("FLIPDIST_MIN_K_LINEAR_PREFIX");
    if (env && *env) {
        int value = std::atoi(env);
        if (value < 0) value = 0;
        if (value > 16) value = 16;
        env_cached = value;
        return env_cached;
    }
    const int n = static_cast<int>(tree.original_nodes.size());
    if (n >= 25) {
        const EmptySShapeStats source = emptySShapeStats(tree);
        const EmptySShapeStats dest = emptySShapeStats(target);
        if (dest.height >= source.height + 2 &&
            dest.avg_depth >= source.avg_depth + 0.6 &&
            dest.branching >= source.branching + 1) {
            return 0;
        }
        return 3;
    }
    return n >= 23 ? 3 : 0;
}

int minKLinearAbortMs(const VectorRangeTreeMap& tree) {
    static int env_cached = INT_MIN;
    if (env_cached != INT_MIN) {
        return env_cached;
    }
    const char* env = std::getenv("FLIPDIST_MIN_K_LINEAR_ABORT_MS");
    if (env && *env) {
        int value = std::atoi(env);
        if (value < 0) value = 0;
        if (value > 10000) value = 10000;
        env_cached = value;
        return env_cached;
    }
    const int n = static_cast<int>(tree.original_nodes.size());
    return n >= 25 ? 150 : 0;
}

void updatePairIncumbentUpper(const Key128& pair_key, int k) {
    if (!incumbentPruneEnabled()) return;
    auto it = g_pair_incumbent_upper.find(pair_key);
    if (it == g_pair_incumbent_upper.end() || k < it->second) {
        g_pair_incumbent_upper[pair_key] = k;
    }
}

int pairIncumbentUpper(const Key128& pair_key) {
    if (!incumbentPruneEnabled()) return std::numeric_limits<int>::max();
    auto it = g_pair_incumbent_upper.find(pair_key);
    if (it == g_pair_incumbent_upper.end()) return std::numeric_limits<int>::max();
    return it->second;
}

const std::unordered_set<RP, RPH>& cachedTargetEdgeSet(const VectorRangeTreeMap& target) {
    const Key128 target_key = makeKeyPair(target, target);
    auto target_it = g_target_edge_set_cache.find(target_key);
    if (target_it == g_target_edge_set_cache.end()) {
        target_it = g_target_edge_set_cache.emplace(target_key, buildTargetSet(target)).first;
    }
    return target_it->second;
}

std::pair<bool, std::pair<int, int>> findFreeEdgeWithCachedTarget(
    const VectorRangeTreeMap& cur,
    const VectorRangeTreeMap& target) {
    const auto& target_edges = cachedTargetEdgeSet(target);
    if (target_edges.empty()) {
        return {false, {-1, -1}};
    }

    std::array<std::pair<int, int>, 64> current_child_ranges{};
    int current_child_range_count = 0;
    bool current_child_ranges_overflow = false;
    bool current_child_ranges_ready = false;
    const auto ensure_current_child_ranges = [&]() {
        if (current_child_ranges_ready) {
            return;
        }
        current_child_ranges_ready = true;
        for (int v : cur.original_nodes) {
            if (!cur.isOriginal(v)) continue;
            const int left = cur.getLeftChild(v);
            const int right = cur.getRightChild(v);
            if (left >= 0 && cur.isOriginal(left)) {
                if (current_child_range_count < static_cast<int>(current_child_ranges.size())) {
                    current_child_ranges[current_child_range_count++] = cur.getRange(left);
                } else {
                    current_child_ranges_overflow = true;
                }
            }
            if (right >= 0 && cur.isOriginal(right)) {
                if (current_child_range_count < static_cast<int>(current_child_ranges.size())) {
                    current_child_ranges[current_child_range_count++] = cur.getRange(right);
                } else {
                    current_child_ranges_overflow = true;
                }
            }
        }
    };
    const auto has_current_child_range = [&](const RP& edge) {
        ensure_current_child_ranges();
        if (current_child_ranges_overflow) {
            return hasEdgeByRange(cur, edge);
        }
        for (int i = 0; i < current_child_range_count; ++i) {
            const auto& range = current_child_ranges[i];
            if (range.first == edge.cs && range.second == edge.ce) {
                return true;
            }
        }
        return false;
    };

    for (int v : cur.original_nodes) {
        if (!cur.isOriginal(v)) continue;

        const int right_child = cur.getRightChild(v);
        if (right_child >= 0 && cur.isOriginal(right_child)) {
            const auto parent_range = cur.getRange(v);
            const auto child_range = cur.getRange(right_child);
            const int left_grandchild = cur.getLeftChild(right_child);
            const auto new_child_range =
                (left_grandchild >= 0 && cur.isOriginal(left_grandchild))
                    ? cur.getRange(left_grandchild)
                    : std::pair<int, int>{parent_range.first, child_range.first};
            const RP edge{parent_range.first,
                          parent_range.second,
                          new_child_range.first,
                          new_child_range.second};
            if (target_edges.count(edge) && !has_current_child_range(edge)) {
                return {true, {v, right_child}};
            }
        }

        const int left_child = cur.getLeftChild(v);
        if (left_child >= 0 && cur.isOriginal(left_child)) {
            const auto parent_range = cur.getRange(v);
            const auto child_range = cur.getRange(left_child);
            const int right_grandchild = cur.getRightChild(left_child);
            const auto new_child_range =
                (right_grandchild >= 0 && cur.isOriginal(right_grandchild))
                    ? cur.getRange(right_grandchild)
                    : std::pair<int, int>{child_range.second, parent_range.second};
            const RP edge{parent_range.first,
                          parent_range.second,
                          new_child_range.first,
                          new_child_range.second};
            if (target_edges.count(edge) && !has_current_child_range(edge)) {
                return {true, {v, left_child}};
            }
        }
    }

    return {false, {-1, -1}};
}

std::uint64_t edgeScoreKey(const std::pair<int, int>& edge) {
    int a = std::min(edge.first, edge.second);
    int b = std::max(edge.first, edge.second);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) ^
           static_cast<std::uint32_t>(b);
}

int collectOutgoingRangeEdgesForNodes(const VectorRangeTreeMap& tree,
                                      const int first_node,
                                      const int second_node,
                                      std::array<RP, 4>& out) {
    int count = 0;
    const auto add_node = [&](const int node) {
        if (node < 0 || !tree.isOriginal(node)) return;
        auto pr = tree.getRange(node);
        for (int child : {tree.getLeftChild(node), tree.getRightChild(node)}) {
            if (child >= 0 && tree.isOriginal(child) && count < static_cast<int>(out.size())) {
                auto cr = tree.getRange(child);
                out[count++] = RP{pr.first, pr.second, cr.first, cr.second};
            }
        }
    };
    add_node(first_node);
    if (second_node != first_node) {
        add_node(second_node);
    }
    return count;
}

int countLocalConflicts(const std::array<RP, 4>& edges,
                        const int edge_count,
                        const std::unordered_set<RP, RPH>& target_edges) {
    int conflicts = 0;
    for (int i = 0; i < edge_count; ++i) {
        if (!target_edges.count(edges[i])) {
            ++conflicts;
        }
    }
    return conflicts;
}

int rotatedConflictCountFromLocalDelta(const int current_conflicts,
                                       const std::array<RP, 4>& old_edges,
                                       const int old_edge_count,
                                       const VectorRangeTreeMap& rotated_tree,
                                       const int first_node,
                                       const int second_node,
                                       const std::unordered_set<RP, RPH>& target_edges) {
    std::array<RP, 4> new_edges{};
    const int new_edge_count =
        collectOutgoingRangeEdgesForNodes(rotated_tree, first_node, second_node, new_edges);
    return current_conflicts -
           countLocalConflicts(old_edges, old_edge_count, target_edges) +
           countLocalConflicts(new_edges, new_edge_count, target_edges);
}

int cachedConflictCountWithKey(const VectorRangeTreeMap& start,
                               const VectorRangeTreeMap& target,
                               const Key128& key) {
    auto it = g_conflict_cache.find(key);
    if (it != g_conflict_cache.end()) {
        if (g_profile.enabled) g_profile.conflict_cache_hits++;
        return it->second;
    }
    int value = 0;
    if (!start.original_inorder.empty() && !target.original_inorder.empty()) {
        const auto& target_edges = cachedTargetEdgeSet(target);
        for (int v : start.original_inorder) {
            if (!start.isOriginal(v)) continue;
            auto pr = start.getRange(v);
            for (int c : {start.getLeftChild(v), start.getRightChild(v)}) {
                if (c >= 0 && start.isOriginal(c)) {
                    auto cr = start.getRange(c);
                    RP edge{pr.first, pr.second, cr.first, cr.second};
                    if (!target_edges.count(edge)) {
                        ++value;
                    }
                }
            }
        }
    }
    g_conflict_cache.emplace(key, value);
    return value;
}

int cachedConflictCount(const VectorRangeTreeMap& start, const VectorRangeTreeMap& target) {
    return cachedConflictCountWithKey(start, target, makeKeyPair(start, target));
}

std::pair<bool, std::pair<int, int>> cachedFindFreeEdge(
    const VectorRangeTreeMap& start, const VectorRangeTreeMap& target, const Key128& key) {
    auto it = g_free_edge_cache.find(key);
    if (it != g_free_edge_cache.end()) {
        if (g_profile.enabled) g_profile.free_edge_cache_hits++;
        return {it->second.has_edge, it->second.edge};
    }
    auto value = findFreeEdgeWithCachedTarget(start, target);
    g_free_edge_cache.emplace(key, CachedFreeEdge{value.first, value.second});
    return value;
}

std::pair<bool, std::pair<int, int>> cachedFindFreeEdge(
    const VectorRangeTreeMap& start, const VectorRangeTreeMap& target) {
    return cachedFindFreeEdge(start, target, makeKeyPair(start, target));
}

bool applyRotationOnEdge(VectorRangeTreeMap& tree, const std::pair<int, int>& edge) {
    auto oriented = orientEdge(tree, edge);
    int parent = oriented.first;
    int child = oriented.second;
    if (parent < 0 || child < 0) {
        return false;
    }
    if (tree.getLeftChild(parent) == child) {
        tree.rotateRight(parent);
        return true;
    }
    if (tree.getRightChild(parent) == child) {
        tree.rotateLeft(parent);
        return true;
    }
    return false;
}

int exactDistanceSmallSubproblem(const VectorRangeTreeMap& start,
                                 const VectorRangeTreeMap& target) {
    if (TreesEqual(start, target)) {
        return 0;
    }
    std::queue<std::pair<VectorRangeTreeMap, int>> q;
    std::unordered_set<Key128, Key128Hash> visited;
    VectorRangeTreeMap start_copy = safeCopyTree(start);
    if (start_copy.original_nodes.empty()) {
        return std::numeric_limits<int>::max() / 4;
    }
    q.push({start_copy, 0});
    visited.insert(makeKeyPair(start_copy, start_copy));

    while (!q.empty()) {
        VectorRangeTreeMap cur = std::move(q.front().first);
        int dist = q.front().second;
        q.pop();
        auto edges = getInternalEdges(cur);
        for (const auto& edge : edges) {
            VectorRangeTreeMap next = safeCopyTree(cur);
            if (next.original_nodes.empty()) continue;
            if (!applyRotationOnEdge(next, edge)) continue;
            Key128 next_key = makeKeyPair(next, next);
            if (!visited.insert(next_key).second) continue;
            if (TreesEqual(next, target)) {
                return dist + 1;
            }
            q.push({std::move(next), dist + 1});
        }
    }

    return std::numeric_limits<int>::max() / 4;
}

int localNoPathLowerBound(const VectorRangeTreeMap& start,
                          const VectorRangeTreeMap& target,
                          int depth) {
    if (depth <= 0 || TreesEqual(start, target)) {
        return 0;
    }
    std::queue<std::pair<VectorRangeTreeMap, int>> q;
    std::unordered_set<Key128, Key128Hash> visited;
    VectorRangeTreeMap start_copy = safeCopyTree(start);
    if (start_copy.original_nodes.empty()) {
        return 0;
    }
    q.push({start_copy, 0});
    visited.insert(makeKeyPair(start_copy, start_copy));

    while (!q.empty()) {
        VectorRangeTreeMap cur = std::move(q.front().first);
        int dist = q.front().second;
        q.pop();
        if (dist >= depth) {
            continue;
        }
        auto edges = getInternalEdges(cur);
        for (const auto& edge : edges) {
            VectorRangeTreeMap next = safeCopyTree(cur);
            if (next.original_nodes.empty()) continue;
            if (!applyRotationOnEdge(next, edge)) continue;
            Key128 next_key = makeKeyPair(next, next);
            if (!visited.insert(next_key).second) continue;
            if (TreesEqual(next, target)) {
                return 0;
            }
            q.push({std::move(next), dist + 1});
        }
    }
    return depth + 1;
}

Key128 makeHeuristicStateKey(
    const VectorRangeTreeMap& start,
    const VectorRangeTreeMap& target,
    HeuristicStateKind kind,
    const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& S_state,
    const std::vector<std::pair<int, int>>& I_state) {
    Key128 key{};
    if (kind == HeuristicStateKind::S) {
        key = makeKeySBase(start, target, S_state);
    } else if (kind == HeuristicStateKind::I) {
        key = makeKeyIBase(start, target, I_state);
    } else {
        key = makeKeyPair(start, target);
    }

    key = mixHeuristicTag(key, 0x11ULL + static_cast<std::uint64_t>(kind));
    key = mixHeuristicTag(key, heuristicConfigTag());
    return key;
}

int heuristicLowerBound(
    const VectorRangeTreeMap& start,
    const VectorRangeTreeMap& target,
    HeuristicStateKind kind,
    const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& S_state,
    const std::vector<std::pair<int, int>>& I_state,
    int conflict_lb_hint) {
    if (!combinedLowerBoundEnabled()) {
        return (conflict_lb_hint >= 0) ? conflict_lb_hint : cachedConflictCount(start, target);
    }
    if (g_profile.enabled) g_profile.heuristic_lb_calls++;
    const bool empty_s_state = (kind == HeuristicStateKind::S && S_state.empty());
    const bool default_bound_config = defaultEmptySBoundConfig();
    if (default_bound_config && empty_s_state) {
        const EmptySBaseKeys empty_s_keys = makeEmptySBaseKeys(start, target);
        const int conflict_lb = (conflict_lb_hint >= 0) ? conflict_lb_hint : cachedConflictCount(start, target);
        return emptySLowerBoundFromKeys(conflict_lb, empty_s_keys.pair_key, empty_s_keys.bounds_key);
    }
    const EmptySKeys empty_s_keys = empty_s_state ? makeEmptySKeys(start, target, 0) : EmptySKeys{};
    Key128 cache_key = empty_s_state
        ? mixHeuristicTag(empty_s_keys.bounds_key, 0x11ULL + static_cast<std::uint64_t>(kind))
        : makeHeuristicStateKey(start, target, kind, S_state, I_state);
    if (empty_s_state) {
        cache_key = mixHeuristicTag(cache_key, heuristicConfigTag());
    }
    auto hit = g_heuristic_bundle_cache.find(cache_key);
    if (hit != g_heuristic_bundle_cache.end()) {
        if (g_profile.enabled) g_profile.heuristic_lb_cache_hits++;
        return hit->second.combined_lb;
    }

    const Key128 pair_key = empty_s_state ? empty_s_keys.pair_key : makeKeyPair(start, target);
    HeuristicBundle bundle;
    bundle.conflict_lb = (conflict_lb_hint >= 0) ? conflict_lb_hint : cachedConflictCount(start, target);
    bundle.state_bounds_lb = requiredBudgetFromBounds(g_kbounds, pair_key);

    if (kind == HeuristicStateKind::S) {
        const Key128 state_key = empty_s_state ? empty_s_keys.bounds_key : makeKeySBase(start, target, S_state);
        bundle.state_bounds_lb = std::max(bundle.state_bounds_lb,
                                          requiredBudgetFromBounds(g_tds_bounds, state_key));
    } else if (kind == HeuristicStateKind::I) {
        const Key128 state_key = makeKeyIBase(start, target, I_state);
        bundle.state_bounds_lb = std::max(bundle.state_bounds_lb,
                                          requiredBudgetFromBounds(g_tdi_bounds, state_key));
    }

    const int small_m = smallExactBoundM();
    if (small_m > 0 &&
        countInternalEdges(start) <= small_m &&
        countInternalEdges(target) <= small_m) {
        if (g_profile.enabled) g_profile.small_exact_lookups++;
        auto it = g_small_exact_dist_cache.find(pair_key);
        if (it != g_small_exact_dist_cache.end()) {
            if (g_profile.enabled) g_profile.small_exact_hits++;
            bundle.small_exact_lb = it->second;
        } else {
            int dist = exactDistanceSmallSubproblem(start, target);
            if (dist >= std::numeric_limits<int>::max() / 8) {
                dist = 0;
            }
            g_small_exact_dist_cache.emplace(pair_key, dist);
            bundle.small_exact_lb = dist;
        }
    }

    const int no_path_depth = localNoPathDepth();
    if (no_path_depth > 0) {
        if (g_profile.enabled) g_profile.local_nopath_lookups++;
        Key128 no_path_key = mixHeuristicTag(pair_key, 0x4E4F504154485F4CULL);
        no_path_key = mixHeuristicTag(no_path_key, static_cast<std::uint64_t>(no_path_depth));
        auto it = g_local_nopath_lb_cache.find(no_path_key);
        if (it != g_local_nopath_lb_cache.end()) {
            if (g_profile.enabled) g_profile.local_nopath_hits++;
            bundle.local_nopath_lb = it->second;
        } else {
            int lb = localNoPathLowerBound(start, target, no_path_depth);
            g_local_nopath_lb_cache.emplace(no_path_key, lb);
            bundle.local_nopath_lb = lb;
        }
    }

    if (branchyCoreLbEnabled() &&
        kind == HeuristicStateKind::S &&
        S_state.empty()) {
        BranchyCoreSummary branchy_summary =
            analyzeBranchyCoreSummary(start, target, branchyCoreExactM());
        if (branchy_summary.reduction_triggered &&
            branchy_summary.exact_bound_used) {
            bundle.branchy_core_lb = branchy_summary.exact_bound;
            if (g_profile.enabled) {
                g_profile.s_empty_branchy_exact_bounds++;
            }
        }
    }

    bundle.combined_lb = std::max({bundle.conflict_lb,
                                   bundle.state_bounds_lb,
                                   bundle.small_exact_lb,
                                   bundle.local_nopath_lb,
                                   bundle.branchy_core_lb});
    g_heuristic_bundle_cache.emplace(cache_key, bundle);
    return bundle.combined_lb;
}

int combinedLowerBoundS(const VectorRangeTreeMap& start,
                        const VectorRangeTreeMap& target,
                        const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& S_state,
                        int conflict_lb) {
    static const std::vector<std::pair<int, int>> empty_i;
    return heuristicLowerBound(start, target, HeuristicStateKind::S, S_state, empty_i, conflict_lb);
}

int combinedLowerBoundI(const VectorRangeTreeMap& start,
                        const VectorRangeTreeMap& target,
                        const std::vector<std::pair<int, int>>& I_state,
                        int conflict_lb) {
    static const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> empty_s;
    return heuristicLowerBound(start, target, HeuristicStateKind::I, empty_s, I_state, conflict_lb);
}

bool evaluateRotationHeuristic(const VectorRangeTreeMap& start,
                               const VectorRangeTreeMap& target,
                               const std::pair<int, int>& edge,
                               int current_conflicts,
                               int& gain_out,
                               int& next_conflicts_out,
                               int& free_hint_out,
                               int& heuristic_lb_out) {
    auto oriented = orientEdge(start, edge);
    int parent = oriented.first;
    int child = oriented.second;
    if (parent < 0 || child < 0) {
        return false;
    }
    VectorRangeTreeMap probe = safeCopyTree(start);
    if (probe.original_nodes.empty()) {
        return false;
    }
    if (probe.getLeftChild(parent) == child) {
        probe.rotateRight(parent);
    } else if (probe.getRightChild(parent) == child) {
        probe.rotateLeft(parent);
    } else {
        return false;
    }

    next_conflicts_out = cachedConflictCount(probe, target);
    gain_out = current_conflicts - next_conflicts_out;
    auto [has_free, _] = cachedFindFreeEdge(probe, target);
    free_hint_out = has_free ? 1 : 0;
    static const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> empty_s;
    heuristic_lb_out = combinedLowerBoundS(probe, target, empty_s, next_conflicts_out);
    return true;
}

void scoreAndOrderSPairs(
    const VectorRangeTreeMap& start,
    const VectorRangeTreeMap& target,
    int current_conflicts,
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& S_pairs,
    std::unordered_map<std::uint64_t, int>& edge_gain,
    std::unordered_map<std::uint64_t, int>& edge_free_hint,
    std::unordered_map<std::uint64_t, int>& edge_lb_hint) {

    struct PairScore {
        std::pair<std::pair<int, int>, std::pair<int, int>> pair;
        int valid_options = 0;
        int best_fscore = std::numeric_limits<int>::max();
        int best_gain = std::numeric_limits<int>::min() / 2;
        int best_conflicts = std::numeric_limits<int>::max();
        int best_free_hint = 0;
    };

    auto scoreEdge = [&](const std::pair<int, int>& e) {
        const std::uint64_t k = edgeScoreKey(e);
        if (edge_gain.find(k) != edge_gain.end()) {
            return;
        }
        int gain = std::numeric_limits<int>::min() / 2;
        int next_conflicts = std::numeric_limits<int>::max();
        int free_hint = 0;
        int child_lb = std::numeric_limits<int>::max() / 2;
        if (evaluateRotationHeuristic(start, target, e, current_conflicts, gain,
                                      next_conflicts, free_hint, child_lb)) {
            edge_gain[k] = gain;
            edge_free_hint[k] = free_hint;
            edge_lb_hint[k] = child_lb;
        } else {
            edge_gain[k] = std::numeric_limits<int>::min() / 2;
            edge_free_hint[k] = 0;
            edge_lb_hint[k] = std::numeric_limits<int>::max() / 2;
        }
    };

    std::vector<PairScore> scored;
    scored.reserve(S_pairs.size());
    for (const auto& p : S_pairs) {
        scoreEdge(p.first);
        scoreEdge(p.second);

        PairScore s;
        s.pair = p;
        const std::uint64_t k1 = edgeScoreKey(p.first);
        const std::uint64_t k2 = edgeScoreKey(p.second);
        const int g1 = edge_gain[k1];
        const int g2 = edge_gain[k2];
        const int lb1 = edge_lb_hint[k1];
        const int lb2 = edge_lb_hint[k2];
        const bool e1_valid = (lb1 < std::numeric_limits<int>::max() / 4);
        const bool e2_valid = (lb2 < std::numeric_limits<int>::max() / 4);
        if (e1_valid) {
            s.valid_options++;
        }
        if (!(p.first == p.second) && e2_valid) {
            s.valid_options++;
        }
        s.best_gain = std::max(g1, g2);
        s.best_free_hint = std::max(edge_free_hint[k1], edge_free_hint[k2]);
        s.best_fscore = 1 + std::min(lb1, lb2);

        int nc = std::numeric_limits<int>::max();
        if (g1 > std::numeric_limits<int>::min() / 4) {
            nc = std::min(nc, current_conflicts - g1);
        }
        if (g2 > std::numeric_limits<int>::min() / 4) {
            nc = std::min(nc, current_conflicts - g2);
        }
        s.best_conflicts = nc;
        scored.push_back(s);
    }

    std::sort(scored.begin(), scored.end(), [](const PairScore& a, const PairScore& b) {
        if (a.valid_options != b.valid_options) return a.valid_options < b.valid_options;
        if (a.best_fscore != b.best_fscore) return a.best_fscore < b.best_fscore;
        if (a.best_gain != b.best_gain) return a.best_gain > b.best_gain;
        if (a.best_free_hint != b.best_free_hint) return a.best_free_hint > b.best_free_hint;
        if (a.best_conflicts != b.best_conflicts) return a.best_conflicts < b.best_conflicts;
        if (a.pair.first != b.pair.first) return a.pair.first < b.pair.first;
        return a.pair.second < b.pair.second;
    });

    S_pairs.clear();
    S_pairs.reserve(scored.size());
    for (const auto& s : scored) {
        S_pairs.push_back(s.pair);
    }
}

bool pairHeuristicOrderingEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env_astar = std::getenv("FLIPDIST_ASTAR_ORDER");
    const char* env_legacy = std::getenv("FLIPDIST_HEURISTIC_PAIR_ORDER");
    if ((env_astar && std::string(env_astar) == "1") ||
        (env_legacy && std::string(env_legacy) == "1")) {
        cached = 1;
    } else {
        cached = 0;
    }
    return cached == 1;
}

bool progressiveSPrefixEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_PROGRESSIVE_S");
    if (!env || std::string(env).empty()) {
        cached = 1;
        return true;
    }
    cached = (std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

int progressiveSPrefixStart() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_PROGRESSIVE_S_START");
    if (!env || std::string(env).empty()) {
        cached = 12;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 1) value = 1;
    if (value > 256) value = 256;
    cached = value;
    return cached;
}

bool emptySCoreReduceEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_EMPTY_S_CORE_REDUCE");
    cached = (env && std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

bool plateauTreeDumpEnabled() {
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const char* env = std::getenv("FLIPDIST_PROFILE_PLATEAU_DUMP_TREES");
    cached = (env && std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

int plateauTreeDumpMinEdges() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) return cached;
    const char* env = std::getenv("FLIPDIST_PROFILE_PLATEAU_DUMP_MIN_EDGES");
    cached = (env && *env) ? std::atoi(env) : 11;
    if (cached < 0) cached = 0;
    return cached;
}

int plateauTreeDumpMaxEdges() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) return cached;
    const char* env = std::getenv("FLIPDIST_PROFILE_PLATEAU_DUMP_MAX_EDGES");
    cached = (env && *env) ? std::atoi(env) : 14;
    if (cached < 0) cached = 0;
    return cached;
}

int plateauTreeDumpLimit() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) return cached;
    const char* env = std::getenv("FLIPDIST_PROFILE_PLATEAU_DUMP_LIMIT");
    cached = (env && *env) ? std::atoi(env) : 12;
    if (cached < 1) cached = 1;
    return cached;
}

std::string plateauTreeDumpBranchPairs() {
    const char* env = std::getenv("FLIPDIST_PROFILE_PLATEAU_DUMP_BRANCH_PAIRS");
    return env ? std::string(env) : std::string();
}

std::string plateauTreeDumpFile() {
    const char* env = std::getenv("FLIPDIST_PROFILE_PLATEAU_DUMP_FILE");
    return env ? std::string(env) : std::string();
}

bool plateauTreeDumpBranchPairMatches(int a, int b) {
    std::string spec = plateauTreeDumpBranchPairs();
    if (spec.empty()) return true;
    const int lo = std::min(a, b);
    const int hi = std::max(a, b);
    std::stringstream ss(spec);
    std::string token;
    while (std::getline(ss, token, ';')) {
        std::stringstream ts(token);
        std::string lhs;
        std::string rhs;
        if (!std::getline(ts, lhs, ',')) continue;
        if (!std::getline(ts, rhs, ',')) continue;
        int x = std::atoi(lhs.c_str());
        int y = std::atoi(rhs.c_str());
        if (std::make_pair(lo, hi) == std::make_pair(std::min(x, y), std::max(x, y))) {
            return true;
        }
    }
    return false;
}

std::string currentTreeToCanonicalString(const VectorRangeTreeMap& tree) {
    std::vector<int> preorder;
    std::vector<int> inorder;
    std::function<void(int)> dfs_pre = [&](int node) {
        if (node < 0 || !tree.isOriginal(node)) return;
        preorder.push_back(node);
        dfs_pre(tree.getLeftChild(node));
        dfs_pre(tree.getRightChild(node));
    };
    std::function<void(int)> dfs_in = [&](int node) {
        if (node < 0 || !tree.isOriginal(node)) return;
        dfs_in(tree.getLeftChild(node));
        inorder.push_back(node);
        dfs_in(tree.getRightChild(node));
    };
    dfs_pre(tree.root);
    dfs_in(tree.root);
    auto join = [](const std::vector<int>& values) {
        std::ostringstream oss;
        for (size_t i = 0; i < values.size(); ++i) {
            if (i) oss << ",";
            oss << values[i];
        }
        return oss.str();
    };
    std::ostringstream out;
    out << "P:" << join(preorder) << ";I:" << join(inorder);
    return out.str();
}

void maybeDumpPlateauTreeState(const Key128& state_key,
                               const VectorRangeTreeMap& start,
                               const VectorRangeTreeMap& target,
                               const TreeStructureProfile& start_profile,
                               const TreeStructureProfile& target_profile,
                               const std::string& motif,
                               int reduced_edges,
                               int conflicts,
                               int k,
                               int legal_children,
                               int plateau_buckets) {
    if (!plateauTreeDumpEnabled()) return;
    if (g_plateau_tree_dump_count >= plateauTreeDumpLimit()) return;
    if (reduced_edges < plateauTreeDumpMinEdges() || reduced_edges > plateauTreeDumpMaxEdges()) return;
    if (!(motif == "branchy" || motif == "branchy_to_broom")) return;
    if (!plateauTreeDumpBranchPairMatches(start_profile.branching_nodes, target_profile.branching_nodes)) return;
    if (!g_plateau_tree_dumped_states.insert(state_key).second) return;

    std::string path = plateauTreeDumpFile();
    if (path.empty()) {
        path = "results/plateau_tree_dumps.jsonl";
    }
    std::ofstream out(path, std::ios::app);
    if (!out) return;
    out << "{"
        << "\"state_hi\":\"" << state_key.hi << "\""
        << ",\"state_lo\":\"" << state_key.lo << "\""
        << ",\"motif\":\"" << motif << "\""
        << ",\"reduced_edges\":" << reduced_edges
        << ",\"conflicts\":" << conflicts
        << ",\"k\":" << k
        << ",\"legal_children\":" << legal_children
        << ",\"plateau_buckets\":" << plateau_buckets
        << ",\"start_branching_nodes\":" << start_profile.branching_nodes
        << ",\"target_branching_nodes\":" << target_profile.branching_nodes
        << ",\"start_max_branch_subtree_edges\":" << start_profile.max_branch_subtree_edges
        << ",\"target_max_branch_subtree_edges\":" << target_profile.max_branch_subtree_edges
        << ",\"start_degree_histogram\":\"" << start_profile.degree_histogram << "\""
        << ",\"target_degree_histogram\":\"" << target_profile.degree_histogram << "\""
        << ",\"start_tree\":\"" << currentTreeToCanonicalString(start) << "\""
        << ",\"target_tree\":\"" << currentTreeToCanonicalString(target) << "\""
        << "}\n";
    ++g_plateau_tree_dump_count;
}

int emptySCoreDbM() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_EMPTY_S_CORE_DB_M");
    if (!env || std::string(env).empty()) {
        cached = 8;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 0) value = 0;
    if (value > 16) value = 16;
    cached = value;
    return cached;
}

bool emptySMotifReduceEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_EMPTY_S_MOTIF_REDUCE");
    cached = (env && std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

int emptySExactMotifM() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_EMPTY_S_EXACT_MOTIF_M");
    if (!env || std::string(env).empty()) {
        cached = 7;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 0) value = 0;
    if (value > 12) value = 12;
    cached = value;
    return cached;
}

bool branchyCoreLbEnabled() {
    static int cached = -1;
    if (cached >= 0) {
        return cached == 1;
    }
    const char* env = std::getenv("FLIPDIST_BRANCHY_CORE_LB");
    cached = (env && std::string(env) == "1") ? 1 : 0;
    return cached == 1;
}

int branchyCoreExactM() {
    static int cached = INT_MIN;
    if (cached != INT_MIN) {
        return cached;
    }
    const char* env = std::getenv("FLIPDIST_BRANCHY_CORE_EXACT_M");
    if (!env || std::string(env).empty()) {
        cached = 8;
        return cached;
    }
    int value = std::atoi(env);
    if (value < 0) value = 0;
    if (value > 12) value = 12;
    cached = value;
    return cached;
}

std::string classifyTreeMotif(const VectorRangeTreeMap& tree) {
    int branching_nodes = 0;
    for (int node : tree.original_nodes) {
        if (!tree.isOriginal(node)) continue;
        int child_count = 0;
        if (tree.getLeftChild(node) >= 0 && tree.isOriginal(tree.getLeftChild(node))) child_count++;
        if (tree.getRightChild(node) >= 0 && tree.isOriginal(tree.getRightChild(node))) child_count++;
        if (child_count >= 2) {
            branching_nodes++;
        }
    }
    if (branching_nodes == 0) {
        return "chain";
    }
    if (branching_nodes == 1) {
        return "broom";
    }
    return "branchy";
}

std::string classifyEmptySCoreMotif(const VectorRangeTreeMap& start,
                                    const VectorRangeTreeMap& target) {
    const std::string left = classifyTreeMotif(start);
    const std::string right = classifyTreeMotif(target);
    if (left == right) {
        return left;
    }
    return left + "_to_" + right;
}

std::string coreSignatureForPair(const VectorRangeTreeMap& start,
                                 const VectorRangeTreeMap& target,
                                 const std::string& motif) {
    const Key128 pair_key = makeKeyPair(start, target);
    return motif + ":" + std::to_string(countInternalEdges(start)) + ":" +
           std::to_string(pair_key.hi) + ":" + std::to_string(pair_key.lo);
}

bool collectBestCommonEdgeCoreSplit(const VectorRangeTreeMap& start,
                                    const VectorRangeTreeMap& target,
                                    CoreDecompositionCandidate& best_out) {
    if (start.original_nodes.empty() || target.original_nodes.empty()) {
        return false;
    }

    auto target_edges = buildTargetSet(target);
    bool found = false;
    CoreDecompositionCandidate best;

    for (int parent : start.original_nodes) {
        if (!start.isOriginal(parent)) continue;
        auto parent_range = start.getRange(parent);
        for (int child : {start.getLeftChild(parent), start.getRightChild(parent)}) {
            if (child < 0 || !start.isOriginal(child)) continue;
            auto child_range = start.getRange(child);
            RP edge{parent_range.first, parent_range.second, child_range.first, child_range.second};
            if (!target_edges.count(edge)) {
                continue;
            }

            auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(start, parent_range, child_range);
            auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(target, parent_range, child_range);
            if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes) {
                continue;
            }

            const int conf1 = cachedConflictCount(A1, B1);
            const int conf2 = cachedConflictCount(A2, B2);
            const int lb1 = std::max(conf1, requiredBudgetFromBounds(g_kbounds, makeKeyPair(A1, B1)));
            const int lb2 = std::max(conf2, requiredBudgetFromBounds(g_kbounds, makeKeyPair(A2, B2)));
            CoreDecompositionCandidate cand;
            cand.A1 = std::move(A1);
            cand.A2 = std::move(A2);
            cand.B1 = std::move(B1);
            cand.B2 = std::move(B2);
            cand.score = lb1 + lb2;
            cand.peak = std::max(lb1, lb2);
            cand.min_side_edges = std::min(countInternalEdges(cand.A1), countInternalEdges(cand.A2));

            if (!found ||
                cand.score > best.score ||
                (cand.score == best.score && cand.peak > best.peak) ||
                (cand.score == best.score && cand.peak == best.peak &&
                 cand.min_side_edges > best.min_side_edges)) {
                best = std::move(cand);
                found = true;
            }
        }
    }

    if (found) {
        best_out = std::move(best);
    }
    return found;
}

bool tryFreeEdgeCoreReduction(const VectorRangeTreeMap& start,
                              const VectorRangeTreeMap& target,
                              FreeEdgeCoreReduction& out) {
    out = FreeEdgeCoreReduction{};
    auto [has_free, free_edge] = cachedFindFreeEdge(start, target);
    if (!has_free) {
        return false;
    }

    VectorRangeTreeMap rotated = safeCopyTree(start);
    if (rotated.original_nodes.empty()) {
        return false;
    }

    int parent = free_edge.first;
    int child = free_edge.second;
    int u = -1;
    int v = -1;
    if (rotated.getLeftChild(parent) == child) {
        rotated.rotateRight(parent);
        u = child;
        v = parent;
    } else if (rotated.getRightChild(parent) == child) {
        rotated.rotateLeft(parent);
        u = child;
        v = parent;
    } else {
        return false;
    }

    out.valid = true;
    if (TreesEqual(rotated, target)) {
        out.solved = true;
        return true;
    }

    auto child_range = rotated.getRange(v);
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> seeded_pairs;
    appendPartnerPairsFromDiagonals(rotated, {{child_range.first, child_range.second}}, seeded_pairs);
    dedupeSPairs(seeded_pairs);

    auto parent_range = rotated.getRange(u);
    auto [left_start, right_start] =
        VectorRangeTreeMap::partitionAlongEdge(rotated, parent_range, child_range);
    auto [left_target, right_target] =
        VectorRangeTreeMap::partitionAlongEdge(target, parent_range, child_range);
    if (left_start.original_nodes != left_target.original_nodes ||
        right_start.original_nodes != right_target.original_nodes) {
        out.valid = false;
        return false;
    }

    auto [left_pairs, right_pairs] = partitionS(seeded_pairs, left_start, right_start);
    if (!left_pairs.empty() || !right_pairs.empty()) {
        return false;
    }

    out.left_start = std::move(left_start);
    out.left_target = std::move(left_target);
    out.right_start = std::move(right_start);
    out.right_target = std::move(right_target);
    return true;
}

EmptySCoreSummary analyzeEmptySCore(const VectorRangeTreeMap& start,
                                    const VectorRangeTreeMap& target) {
    const Key128 pair_key = makeKeyPair(start, target);
    if (g_profile.enabled) {
        g_profile.s_empty_core_lookups++;
    }

    auto cache_it = g_empty_s_core_cache.find(pair_key);
    if (cache_it != g_empty_s_core_cache.end()) {
        if (g_profile.enabled) {
            g_profile.s_empty_core_cache_hits++;
        }
        return cache_it->second;
    }

    EmptySCoreSummary summary;
    summary.motif = classifyEmptySCoreMotif(start, target);
    summary.signature = coreSignatureForPair(start, target, summary.motif);

    throwIfProfileAbort();

    if (TreesEqual(start, target)) {
        summary.fully_exact = true;
        summary.signature = coreSignatureForPair(start, target, "equal");
        g_empty_s_core_cache.emplace(pair_key, summary);
        return summary;
    }

    CoreDecompositionCandidate split;
    if (collectBestCommonEdgeCoreSplit(start, target, split)) {
        if (g_profile.enabled) {
            g_profile.s_empty_core_common_splits++;
        }
        EmptySCoreSummary left = analyzeEmptySCore(split.A1, split.B1);
        EmptySCoreSummary right = analyzeEmptySCore(split.A2, split.B2);
        summary.lower_bound = left.lower_bound + right.lower_bound;
        summary.common_splits = 1 + left.common_splits + right.common_splits;
        summary.free_reductions = left.free_reductions + right.free_reductions;
        summary.fully_exact = left.fully_exact && right.fully_exact;
        if (summary.fully_exact) {
            summary.exact_total = left.exact_total + right.exact_total;
        }
        summary.motif = "decomposed";
        summary.signature = coreSignatureForPair(start, target, summary.motif);
        g_empty_s_core_cache.emplace(pair_key, summary);
        return summary;
    }

    FreeEdgeCoreReduction free_reduction;
    if (tryFreeEdgeCoreReduction(start, target, free_reduction)) {
        if (g_profile.enabled) {
            g_profile.s_empty_core_free_reductions++;
        }
        summary.motif = "free_reduced";
        summary.signature = coreSignatureForPair(start, target, summary.motif);
        summary.free_reductions = 1;
        if (free_reduction.solved) {
            summary.lower_bound = 1;
            summary.exact_total = 1;
            summary.fully_exact = true;
            g_empty_s_core_cache.emplace(pair_key, summary);
            return summary;
        }

        EmptySCoreSummary left = analyzeEmptySCore(free_reduction.left_start, free_reduction.left_target);
        EmptySCoreSummary right = analyzeEmptySCore(free_reduction.right_start, free_reduction.right_target);
        summary.lower_bound = 1 + left.lower_bound + right.lower_bound;
        summary.common_splits = left.common_splits + right.common_splits;
        summary.free_reductions += left.free_reductions + right.free_reductions;
        summary.fully_exact = left.fully_exact && right.fully_exact;
        if (summary.fully_exact) {
            summary.exact_total = 1 + left.exact_total + right.exact_total;
        }
        g_empty_s_core_cache.emplace(pair_key, summary);
        return summary;
    }

    summary.lower_bound = cachedConflictCount(start, target);
    const int max_edges = emptySCoreDbM();
    const int internal_edges = countInternalEdges(start);
    if (max_edges > 0 && internal_edges <= max_edges) {
        auto exact_it = g_empty_s_core_exact_cache.find(pair_key);
        if (exact_it != g_empty_s_core_exact_cache.end()) {
            if (g_profile.enabled) {
                g_profile.s_empty_core_exact_hits++;
            }
            summary.exact_total = exact_it->second;
            summary.lower_bound = std::max(summary.lower_bound, summary.exact_total);
            summary.fully_exact = true;
        } else {
            int exact_dist = FindRotationDistance(start, target);
            if (g_profile.enabled) {
                g_profile.s_empty_core_exact_hits++;
            }
            g_empty_s_core_exact_cache.emplace(pair_key, exact_dist);
            summary.exact_total = exact_dist;
            summary.lower_bound = std::max(summary.lower_bound, exact_dist);
            summary.fully_exact = true;
        }
    }

    g_empty_s_core_cache.emplace(pair_key, summary);
    return summary;
}

bool branchyResidualEligible(const std::string& motif) {
    return motif == "branchy" || motif == "branchy_to_broom";
}

Key128 branchyCoreSummaryKey(const VectorRangeTreeMap& start,
                             const VectorRangeTreeMap& target,
                             int exact_m) {
    Key128 key = makeKeyPair(start, target);
    key = mixHeuristicTag(key, 0x6272616e63687963ULL);
    key.hi ^= static_cast<std::uint64_t>(exact_m);
    return key;
}

int exactBranchyCoreDistance(const VectorRangeTreeMap& start,
                             const VectorRangeTreeMap& target,
                             bool count_profile) {
    const bool count_hits = count_profile || (g_profile.enabled && branchyCoreLbEnabled());
    const Key128 pair_key = makeKeyPair(start, target);
    auto cache_it = g_branchy_core_exact_cache.find(pair_key);
    if (cache_it != g_branchy_core_exact_cache.end()) {
        if (g_profile.enabled && count_hits) {
            g_profile.s_empty_branchy_exact_cache_hits++;
        }
        return cache_it->second;
    }
    int exact_dist = FindRotationDistance(start, target);
    g_branchy_core_exact_cache.emplace(pair_key, exact_dist);
    return exact_dist;
}

struct BranchyPeelCandidate {
    VectorRangeTreeMap residual_start;
    VectorRangeTreeMap residual_target;
    int peel_cost = 0;
    int residual_edges = 0;
    int residual_conflicts = 0;
    bool from_free = false;
    std::string residual_motif = "unknown";
    std::string residual_signature;
};

bool betterBranchyPeelCandidate(const BranchyPeelCandidate& a,
                                const BranchyPeelCandidate& b) {
    const bool a_branchy = branchyResidualEligible(a.residual_motif);
    const bool b_branchy = branchyResidualEligible(b.residual_motif);
    if (a_branchy != b_branchy) {
        return a_branchy > b_branchy;
    }
    if (a.residual_edges != b.residual_edges) {
        return a.residual_edges < b.residual_edges;
    }
    if (a.peel_cost != b.peel_cost) {
        return a.peel_cost < b.peel_cost;
    }
    if (a.from_free != b.from_free) {
        return a.from_free < b.from_free;
    }
    return a.residual_signature < b.residual_signature;
}

void maybeUpdateBestBranchyCandidate(const VectorRangeTreeMap& peeled_start,
                                     const VectorRangeTreeMap& peeled_target,
                                     const VectorRangeTreeMap& residual_start,
                                     const VectorRangeTreeMap& residual_target,
                                     int base_cost,
                                     bool from_free,
                                     int exact_m,
                                     int current_edges,
                                     bool& found,
                                     BranchyPeelCandidate& best) {
    if (peeled_start.original_nodes != peeled_target.original_nodes ||
        residual_start.original_nodes != residual_target.original_nodes) {
        return;
    }
    const int peeled_edges = countInternalEdges(peeled_start);
    const int residual_edges = countInternalEdges(residual_start);
    if (peeled_edges > exact_m || residual_edges >= current_edges) {
        return;
    }

    int peeled_exact = exactBranchyCoreDistance(peeled_start, peeled_target, false);
    BranchyPeelCandidate cand;
    cand.residual_start = residual_start;
    cand.residual_target = residual_target;
    cand.peel_cost = base_cost + peeled_exact;
    cand.residual_edges = residual_edges;
    cand.residual_conflicts = cachedConflictCount(residual_start, residual_target);
    cand.from_free = from_free;
    cand.residual_motif = classifyEmptySCoreMotif(residual_start, residual_target);
    cand.residual_signature = coreSignatureForPair(residual_start, residual_target, cand.residual_motif);

    if (!(branchyResidualEligible(cand.residual_motif) ||
          TreesEqual(residual_start, residual_target) ||
          residual_edges <= exact_m)) {
        return;
    }

    if (!found || betterBranchyPeelCandidate(cand, best)) {
        best = std::move(cand);
        found = true;
    }
}

bool collectBestBranchyPeelCandidate(const VectorRangeTreeMap& start,
                                     const VectorRangeTreeMap& target,
                                     int exact_m,
                                     BranchyPeelCandidate& best_out) {
    bool found = false;
    BranchyPeelCandidate best;
    const int current_edges = countInternalEdges(start);
    auto target_edges = buildTargetSet(target);

    for (int parent : start.original_nodes) {
        if (!start.isOriginal(parent)) continue;
        auto parent_range = start.getRange(parent);
        for (int child : {start.getLeftChild(parent), start.getRightChild(parent)}) {
            if (child < 0 || !start.isOriginal(child)) continue;
            auto child_range = start.getRange(child);
            RP edge{parent_range.first, parent_range.second, child_range.first, child_range.second};
            if (!target_edges.count(edge)) {
                continue;
            }

            auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(start, parent_range, child_range);
            auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(target, parent_range, child_range);
            if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes) {
                continue;
            }

            maybeUpdateBestBranchyCandidate(A1, B1, A2, B2, 0, false, exact_m, current_edges, found, best);
            maybeUpdateBestBranchyCandidate(A2, B2, A1, B1, 0, false, exact_m, current_edges, found, best);
        }
    }

    FreeEdgeCoreReduction free_reduction;
    if (tryFreeEdgeCoreReduction(start, target, free_reduction)) {
        if (free_reduction.solved) {
            BranchyPeelCandidate cand;
            cand.peel_cost = 1;
            cand.residual_edges = 0;
            cand.residual_conflicts = 0;
            cand.from_free = true;
            cand.residual_motif = "equal";
            cand.residual_signature = coreSignatureForPair(start, target, "equal");
            if (!found || betterBranchyPeelCandidate(cand, best)) {
                best = std::move(cand);
                found = true;
            }
        } else {
            maybeUpdateBestBranchyCandidate(free_reduction.left_start, free_reduction.left_target,
                                            free_reduction.right_start, free_reduction.right_target,
                                            1, true, exact_m, current_edges, found, best);
            maybeUpdateBestBranchyCandidate(free_reduction.right_start, free_reduction.right_target,
                                            free_reduction.left_start, free_reduction.left_target,
                                            1, true, exact_m, current_edges, found, best);
        }
    }

    if (found) {
        best_out = std::move(best);
    }
    return found;
}

BranchyCoreSummary analyzeBranchyCoreSummary(const VectorRangeTreeMap& start,
                                             const VectorRangeTreeMap& target,
                                             int exact_m) {
    const Key128 summary_key = branchyCoreSummaryKey(start, target, exact_m);
    auto cache_it = g_branchy_core_summary_cache.find(summary_key);
    if (cache_it != g_branchy_core_summary_cache.end()) {
        return cache_it->second;
    }

    BranchyCoreSummary summary;
    VectorRangeTreeMap cur_start = safeCopyTree(start);
    VectorRangeTreeMap cur_target = safeCopyTree(target);
    if (cur_start.original_nodes.empty() || cur_target.original_nodes.empty()) {
        summary.residual_edges = 0;
        summary.residual_signature = coreSignatureForPair(start, target, "empty");
        g_branchy_core_summary_cache.emplace(summary_key, summary);
        return summary;
    }

    while (true) {
        throwIfProfileAbort();
        if (TreesEqual(cur_start, cur_target)) {
            summary.residual_edges = 0;
            summary.residual_conflicts = 0;
            summary.residual_motif = "equal";
            summary.residual_signature = coreSignatureForPair(cur_start, cur_target, summary.residual_motif);
            summary.exact_bound = summary.forced_cost;
            summary.exact_bound_used = true;
            g_branchy_core_summary_cache.emplace(summary_key, summary);
            return summary;
        }

        BranchyPeelCandidate best;
        if (!collectBestBranchyPeelCandidate(cur_start, cur_target, exact_m, best)) {
            break;
        }

        summary.reduction_triggered = true;
        summary.forced_cost += best.peel_cost;
        if (best.residual_edges == 0) {
            summary.residual_edges = 0;
            summary.residual_conflicts = 0;
            summary.residual_motif = "equal";
            summary.residual_signature = best.residual_signature;
            summary.exact_bound = summary.forced_cost;
            summary.exact_bound_used = true;
            g_branchy_core_summary_cache.emplace(summary_key, summary);
            return summary;
        }
        cur_start = std::move(best.residual_start);
        cur_target = std::move(best.residual_target);
    }

    summary.residual_edges = countInternalEdges(cur_start);
    summary.residual_conflicts = cachedConflictCount(cur_start, cur_target);
    summary.residual_motif = classifyEmptySCoreMotif(cur_start, cur_target);
    summary.residual_signature = coreSignatureForPair(cur_start, cur_target, summary.residual_motif);
    if (summary.reduction_triggered &&
        summary.residual_edges <= exact_m &&
        branchyResidualEligible(summary.residual_motif)) {
        int residual_exact = exactBranchyCoreDistance(cur_start, cur_target, false);
        summary.exact_bound = summary.forced_cost + residual_exact;
        summary.exact_bound_used = true;
    }

    g_branchy_core_summary_cache.emplace(summary_key, summary);
    return summary;
}

bool motifPeelEligible(const std::string& motif) {
    return motif == "chain" || motif == "broom" ||
           motif == "chain_to_broom" || motif == "broom_to_chain";
}

bool motifExactFallbackEligible(const std::string& motif) {
    return motif == "chain" || motif == "broom" ||
           motif == "broom_to_chain" || motif == "branchy_to_broom";
}

Key128 emptySMotifSummaryKey(const VectorRangeTreeMap& start,
                             const VectorRangeTreeMap& target,
                             int motif_m) {
    Key128 key = makeKeyPair(start, target);
    key = mixHeuristicTag(key, 0x6d6f74696673756dULL);
    key.hi ^= static_cast<std::uint64_t>(motif_m);
    return key;
}

TreeStructureProfile describeTreeStructure(const VectorRangeTreeMap& tree) {
    TreeStructureProfile profile;
    if (tree.root < 0 || !tree.isOriginal(tree.root) || tree.original_nodes.empty()) {
        profile.motif = "empty";
        return profile;
    }

    std::array<int, 3> degree_hist{0, 0, 0};

    struct NodeStats {
        int node_count = 0;
        int branching_count = 0;
    };

    std::function<NodeStats(int)> dfs = [&](int node) -> NodeStats {
        if (node < 0 || !tree.isOriginal(node)) {
            return {};
        }

        NodeStats left = dfs(tree.getLeftChild(node));
        NodeStats right = dfs(tree.getRightChild(node));
        int child_count = 0;
        if (left.node_count > 0) child_count++;
        if (right.node_count > 0) child_count++;
        degree_hist[std::min(child_count, 2)]++;

        NodeStats cur;
        cur.node_count = 1 + left.node_count + right.node_count;
        cur.branching_count = left.branching_count + right.branching_count + (child_count >= 2 ? 1 : 0);

        if (child_count >= 2) {
            profile.branching_nodes++;
            if (left.node_count > 0) {
                int left_edges = left.node_count - 1;
                profile.max_branch_subtree_edges =
                    std::max(profile.max_branch_subtree_edges, left_edges);
                if (left.branching_count == 0) {
                    profile.chain_arm_count++;
                } else if (left.branching_count == 1) {
                    profile.broom_arm_count++;
                }
            }
            if (right.node_count > 0) {
                int right_edges = right.node_count - 1;
                profile.max_branch_subtree_edges =
                    std::max(profile.max_branch_subtree_edges, right_edges);
                if (right.branching_count == 0) {
                    profile.chain_arm_count++;
                } else if (right.branching_count == 1) {
                    profile.broom_arm_count++;
                }
            }
        }

        return cur;
    };

    dfs(tree.root);
    std::ostringstream hist;
    hist << "0:" << degree_hist[0] << "|1:" << degree_hist[1] << "|2:" << degree_hist[2];
    profile.degree_histogram = hist.str();
    profile.motif = classifyTreeMotif(tree);
    return profile;
}

int exactEmptySMotifDistance(const VectorRangeTreeMap& start,
                             const VectorRangeTreeMap& target,
                             bool count_profile) {
    const Key128 pair_key = makeKeyPair(start, target);
    auto cache_it = g_empty_s_motif_exact_cache.find(pair_key);
    if (cache_it != g_empty_s_motif_exact_cache.end()) {
        if (g_profile.enabled && count_profile) {
            g_profile.s_empty_motif_exact_cache_hits++;
        }
        return cache_it->second;
    }
    if (g_profile.enabled && count_profile) {
        g_profile.s_empty_motif_exact_fallbacks++;
    }
    int exact_dist = FindRotationDistance(start, target);
    g_empty_s_motif_exact_cache.emplace(pair_key, exact_dist);
    return exact_dist;
}

void maybeAppendMotifPeelLayout(const VectorRangeTreeMap& peeled_start,
                                const VectorRangeTreeMap& peeled_target,
                                const VectorRangeTreeMap& residual_start,
                                const VectorRangeTreeMap& residual_target,
                                int base_cost,
                                const std::string& kind,
                                const std::string& articulation_edge,
                                int motif_m,
                                std::unordered_set<Key128, Key128Hash>& seen_residuals,
                                std::vector<MotifPeelLayout>& out) {
    if (peeled_start.original_nodes != peeled_target.original_nodes ||
        residual_start.original_nodes != residual_target.original_nodes) {
        return;
    }

    int peeled_edges = countInternalEdges(peeled_start);
    if (peeled_edges > motif_m) {
        return;
    }

    const std::string peel_motif = classifyEmptySCoreMotif(peeled_start, peeled_target);
    if (!motifPeelEligible(peel_motif)) {
        return;
    }

    Key128 dedupe_key = makeKeyPair(residual_start, residual_target);
    dedupe_key = mixHeuristicTag(dedupe_key,
                                 kind.find("free") == 0 ? 0x667265656d6f7469ULL
                                                        : 0x636f6d6d6f6e6d6fULL);
    if (!seen_residuals.insert(dedupe_key).second) {
        return;
    }

    MotifPeelLayout layout;
    layout.peeled_start = peeled_start;
    layout.peeled_target = peeled_target;
    layout.residual_start = residual_start;
    layout.residual_target = residual_target;
    layout.base_cost = base_cost;
    layout.peeled_edges = peeled_edges;
    layout.residual_edges = countInternalEdges(residual_start);
    layout.peel_motif = peel_motif;
    layout.kind = kind;
    layout.articulation_edge = articulation_edge;
    out.push_back(std::move(layout));
}

void collectCommonEdgeMotifLayouts(const VectorRangeTreeMap& start,
                                   const VectorRangeTreeMap& target,
                                   int motif_m,
                                   std::unordered_set<Key128, Key128Hash>& seen_residuals,
                                   std::vector<MotifPeelLayout>& out) {
    if (start.original_nodes.empty() || target.original_nodes.empty()) {
        return;
    }

    auto target_edges = buildTargetSet(target);
    for (int parent : start.original_nodes) {
        if (!start.isOriginal(parent)) continue;
        auto parent_range = start.getRange(parent);
        for (int child : {start.getLeftChild(parent), start.getRightChild(parent)}) {
            if (child < 0 || !start.isOriginal(child)) continue;
            auto child_range = start.getRange(child);
            RP edge{parent_range.first, parent_range.second, child_range.first, child_range.second};
            if (!target_edges.count(edge)) {
                continue;
            }

            auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(start, parent_range, child_range);
            auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(target, parent_range, child_range);
            if (A1.original_nodes != B1.original_nodes || A2.original_nodes != B2.original_nodes) {
                continue;
            }

            std::ostringstream edge_label;
            edge_label << parent_range.first << "-" << parent_range.second
                       << "|" << child_range.first << "-" << child_range.second;

            maybeAppendMotifPeelLayout(A1, B1, A2, B2, 0, "common",
                                       edge_label.str(), motif_m, seen_residuals, out);
            maybeAppendMotifPeelLayout(A2, B2, A1, B1, 0, "common",
                                       edge_label.str(), motif_m, seen_residuals, out);
        }
    }
}

bool collectFreeEdgeMotifLayouts(const VectorRangeTreeMap& start,
                                 const VectorRangeTreeMap& target,
                                 int motif_m,
                                 std::unordered_set<Key128, Key128Hash>& seen_residuals,
                                 std::vector<MotifPeelLayout>& out) {
    FreeEdgeCoreReduction reduction;
    if (!tryFreeEdgeCoreReduction(start, target, reduction)) {
        return false;
    }

    if (reduction.solved) {
        MotifPeelLayout layout;
        layout.base_cost = 1;
        layout.solves_entire_state = true;
        layout.kind = "free";
        out.push_back(std::move(layout));
        return true;
    }

    maybeAppendMotifPeelLayout(reduction.left_start, reduction.left_target,
                               reduction.right_start, reduction.right_target,
                               1, "free", "free", motif_m, seen_residuals, out);
    maybeAppendMotifPeelLayout(reduction.right_start, reduction.right_target,
                               reduction.left_start, reduction.left_target,
                               1, "free", "free", motif_m, seen_residuals, out);
    return false;
}

std::vector<MotifPeelLayout> collectEmptySMotifPeelLayouts(const VectorRangeTreeMap& start,
                                                           const VectorRangeTreeMap& target,
                                                           int motif_m,
                                                           bool& solved_by_free_edge) {
    solved_by_free_edge = false;
    std::unordered_set<Key128, Key128Hash> seen_residuals;
    std::vector<MotifPeelLayout> layouts;
    collectCommonEdgeMotifLayouts(start, target, motif_m, seen_residuals, layouts);
    solved_by_free_edge = collectFreeEdgeMotifLayouts(start, target, motif_m, seen_residuals, layouts);
    std::sort(layouts.begin(), layouts.end(), [](const MotifPeelLayout& a, const MotifPeelLayout& b) {
        if (a.solves_entire_state != b.solves_entire_state) {
            return a.solves_entire_state > b.solves_entire_state;
        }
        if (a.residual_edges != b.residual_edges) {
            return a.residual_edges < b.residual_edges;
        }
        if (a.peeled_edges != b.peeled_edges) {
            return a.peeled_edges < b.peeled_edges;
        }
        if (a.base_cost != b.base_cost) {
            return a.base_cost < b.base_cost;
        }
        if (a.kind != b.kind) {
            return a.kind < b.kind;
        }
        return a.peel_motif < b.peel_motif;
    });
    return layouts;
}

bool betterArticulationArmReduction(const ArticulationArmReduction& a,
                                    const ArticulationArmReduction& b) {
    if (a.solved != b.solved) {
        return a.solved > b.solved;
    }
    if (a.residual_edges != b.residual_edges) {
        return a.residual_edges < b.residual_edges;
    }
    if (a.forced_cost != b.forced_cost) {
        return a.forced_cost < b.forced_cost;
    }
    if (a.peeled_edges != b.peeled_edges) {
        return a.peeled_edges > b.peeled_edges;
    }
    return a.residual_motif < b.residual_motif;
}

bool tryArticulationArmReduction(const VectorRangeTreeMap& start,
                                 const VectorRangeTreeMap& target,
                                 int k,
                                 ArticulationArmReduction& out) {
    out = ArticulationArmReduction{};
    const int arm_m = articulationArmMaxEdges();
    if (arm_m <= 0 || start.original_nodes.empty() || target.original_nodes.empty()) {
        return false;
    }

    std::unordered_set<Key128, Key128Hash> seen_residuals;
    std::vector<MotifPeelLayout> layouts;
    collectCommonEdgeMotifLayouts(start, target, arm_m, seen_residuals, layouts);
    if (layouts.empty()) {
        return false;
    }

    const int current_edges = countInternalEdges(start);
    const int current_conflicts = cachedConflictCount(start, target);
    const std::string current_motif = classifyEmptySCoreMotif(start, target);
    const std::string current_signature = coreSignatureForPair(start, target, current_motif);
    bool found = false;
    ArticulationArmReduction best;

    auto logCandidate = [&](const MotifPeelLayout& layout,
                            const std::string& reject_reason,
                            int exact_side_cost,
                            int forced_cost,
                            int residual_conflicts,
                            const std::string& residual_motif,
                            const TreeStructureProfile* peeled_start_profile,
                            const TreeStructureProfile* peeled_target_profile,
                            const TreeStructureProfile* residual_start_profile,
                            const TreeStructureProfile* residual_target_profile,
                            bool selected) {
        if (!articulationArmLogEnabled()) {
            return;
        }
        std::ostringstream json;
        json << "{";
        json << "\"state_signature\":\"" << jsonEscape(current_signature) << "\"";
        json << ",\"state_motif\":\"" << jsonEscape(current_motif) << "\"";
        json << ",\"current_edges\":" << current_edges;
        json << ",\"current_conflicts\":" << current_conflicts;
        json << ",\"k\":" << k;
        json << ",\"kind\":\"" << jsonEscape(layout.kind) << "\"";
        json << ",\"articulation_edge\":\"" << jsonEscape(layout.articulation_edge) << "\"";
        json << ",\"peeled_edges\":" << layout.peeled_edges;
        json << ",\"residual_edges\":" << layout.residual_edges;
        json << ",\"exact_side_cost\":" << exact_side_cost;
        json << ",\"forced_cost\":" << forced_cost;
        json << ",\"residual_conflicts\":" << residual_conflicts;
        json << ",\"peeled_motif\":\"" << jsonEscape(layout.peel_motif) << "\"";
        json << ",\"residual_motif\":\"" << jsonEscape(residual_motif) << "\"";
        json << ",\"reject_reason\":\"" << jsonEscape(reject_reason) << "\"";
        json << ",\"selected\":" << (selected ? "true" : "false");
        if (peeled_start_profile) {
            json << ",\"peeled_start_branching\":" << peeled_start_profile->branching_nodes;
            json << ",\"peeled_start_hist\":\"" << jsonEscape(peeled_start_profile->degree_histogram) << "\"";
            json << ",\"peeled_start_tree_motif\":\"" << jsonEscape(peeled_start_profile->motif) << "\"";
        }
        if (peeled_target_profile) {
            json << ",\"peeled_target_branching\":" << peeled_target_profile->branching_nodes;
            json << ",\"peeled_target_hist\":\"" << jsonEscape(peeled_target_profile->degree_histogram) << "\"";
            json << ",\"peeled_target_tree_motif\":\"" << jsonEscape(peeled_target_profile->motif) << "\"";
        }
        if (residual_start_profile) {
            json << ",\"residual_start_branching\":" << residual_start_profile->branching_nodes;
            json << ",\"residual_start_hist\":\"" << jsonEscape(residual_start_profile->degree_histogram) << "\"";
            json << ",\"residual_start_tree_motif\":\"" << jsonEscape(residual_start_profile->motif) << "\"";
        }
        if (residual_target_profile) {
            json << ",\"residual_target_branching\":" << residual_target_profile->branching_nodes;
            json << ",\"residual_target_hist\":\"" << jsonEscape(residual_target_profile->degree_histogram) << "\"";
            json << ",\"residual_target_tree_motif\":\"" << jsonEscape(residual_target_profile->motif) << "\"";
        }
        json << "}";
        appendArticulationArmLogLine(json.str());
    };

    for (const auto& layout : layouts) {
        int exact_side_cost = -1;
        int forced_cost = -1;
        int residual_conflicts = -1;
        std::string residual_motif = "unknown";
        TreeStructureProfile peeled_start_profile;
        TreeStructureProfile peeled_target_profile;
        TreeStructureProfile residual_start_profile;
        TreeStructureProfile residual_target_profile;
        bool have_peeled_profiles = false;
        bool have_residual_profiles = false;
        auto reject = [&](const std::string& reason) {
            logCandidate(layout,
                         reason,
                         exact_side_cost,
                         forced_cost,
                         residual_conflicts,
                         residual_motif,
                         have_peeled_profiles ? &peeled_start_profile : nullptr,
                         have_peeled_profiles ? &peeled_target_profile : nullptr,
                         have_residual_profiles ? &residual_start_profile : nullptr,
                         have_residual_profiles ? &residual_target_profile : nullptr,
                         false);
        };
        if (layout.solves_entire_state || layout.kind != "common") {
            reject("not_common_or_solves");
            continue;
        }
        if (layout.peeled_edges <= 0 || layout.peeled_edges > arm_m ||
            layout.residual_edges >= current_edges) {
            reject("size_filter");
            continue;
        }
        const std::string peeled_left = classifyTreeMotif(layout.peeled_start);
        const std::string peeled_right = classifyTreeMotif(layout.peeled_target);
        const bool simple_arm =
            (peeled_left == "chain" || peeled_left == "broom") &&
            (peeled_right == "chain" || peeled_right == "broom");
        if (!simple_arm) {
            reject("peeled_not_simple");
            continue;
        }

        peeled_start_profile = describeTreeStructure(layout.peeled_start);
        peeled_target_profile = describeTreeStructure(layout.peeled_target);
        have_peeled_profiles = true;

        residual_motif =
            classifyEmptySCoreMotif(layout.residual_start, layout.residual_target);
        if (!branchyResidualEligible(residual_motif)) {
            reject("residual_not_branchy");
            continue;
        }

        residual_start_profile = describeTreeStructure(layout.residual_start);
        residual_target_profile = describeTreeStructure(layout.residual_target);
        have_residual_profiles = true;

        exact_side_cost = exactEmptySMotifDistance(layout.peeled_start, layout.peeled_target, true);
        forced_cost = layout.base_cost + exact_side_cost;
        if (forced_cost > k) {
            reject("forced_cost_exceeds_k");
            continue;
        }

        residual_conflicts = cachedConflictCount(layout.residual_start, layout.residual_target);

        ArticulationArmReduction cand;
        cand.valid = true;
        cand.solved = TreesEqual(layout.residual_start, layout.residual_target);
        cand.residual_start = layout.residual_start;
        cand.residual_target = layout.residual_target;
        cand.exact_side_cost = exact_side_cost;
        cand.forced_cost = forced_cost;
        cand.peeled_edges = layout.peeled_edges;
        cand.residual_edges = layout.residual_edges;
        cand.residual_motif = residual_motif;
        cand.articulation_edge = layout.articulation_edge;
        cand.peel_motif = layout.peel_motif;

        logCandidate(layout,
                     "candidate",
                     exact_side_cost,
                     forced_cost,
                     residual_conflicts,
                     residual_motif,
                     &peeled_start_profile,
                     &peeled_target_profile,
                     &residual_start_profile,
                     &residual_target_profile,
                     false);
        if (!found || betterArticulationArmReduction(cand, best)) {
            best = std::move(cand);
            found = true;
        }
    }

    if (!found) {
        return false;
    }
    if (articulationArmLogEnabled()) {
        TreeStructureProfile best_residual_start = describeTreeStructure(best.residual_start);
        TreeStructureProfile best_residual_target = describeTreeStructure(best.residual_target);
        MotifPeelLayout picked;
        picked.kind = "common";
        picked.articulation_edge = best.articulation_edge;
        picked.peeled_edges = best.peeled_edges;
        picked.residual_edges = best.residual_edges;
        picked.peel_motif = best.peel_motif;
        logCandidate(picked,
                     "selected",
                     best.exact_side_cost,
                     best.forced_cost,
                     cachedConflictCount(best.residual_start, best.residual_target),
                     best.residual_motif,
                     nullptr,
                     nullptr,
                     &best_residual_start,
                     &best_residual_target,
                     true);
    }
    out = std::move(best);
    return true;
}

EmptySMotifSummary analyzeEmptySMotifSummary(const VectorRangeTreeMap& start,
                                             const VectorRangeTreeMap& target,
                                             int motif_m) {
    const Key128 summary_key = emptySMotifSummaryKey(start, target, motif_m);
    auto cache_it = g_empty_s_motif_summary_cache.find(summary_key);
    if (cache_it != g_empty_s_motif_summary_cache.end()) {
        return cache_it->second;
    }

    EmptySMotifSummary summary;
    summary.start_profile = describeTreeStructure(start);
    summary.target_profile = describeTreeStructure(target);
    summary.chain_arm_count =
        summary.start_profile.chain_arm_count + summary.target_profile.chain_arm_count;
    summary.broom_arm_count =
        summary.start_profile.broom_arm_count + summary.target_profile.broom_arm_count;

    if (TreesEqual(start, target)) {
        summary.reduced_core_internal_edges = 0;
        summary.motif = "equal";
        summary.signature = coreSignatureForPair(start, target, summary.motif);
        g_empty_s_motif_summary_cache.emplace(summary_key, summary);
        return summary;
    }

    summary.reduced_core_internal_edges = countInternalEdges(start);
    summary.motif = classifyEmptySCoreMotif(start, target);
    summary.signature = coreSignatureForPair(start, target, summary.motif);

    bool solved_by_free_edge = false;
    auto layouts = collectEmptySMotifPeelLayouts(start, target, motif_m, solved_by_free_edge);
    if (solved_by_free_edge) {
        summary.reduced_core_internal_edges = 0;
        summary.motif = "free_solved";
        summary.signature = coreSignatureForPair(start, target, summary.motif);
        g_empty_s_motif_summary_cache.emplace(summary_key, summary);
        return summary;
    }

    for (const auto& layout : layouts) {
        if (layout.solves_entire_state) {
            summary.reduced_core_internal_edges = 0;
            summary.motif = "free_solved";
            summary.signature = coreSignatureForPair(start, target, summary.motif);
            break;
        }
        EmptySMotifSummary residual =
            analyzeEmptySMotifSummary(layout.residual_start, layout.residual_target, motif_m);
        if (residual.reduced_core_internal_edges < summary.reduced_core_internal_edges) {
            summary.reduced_core_internal_edges = residual.reduced_core_internal_edges;
            summary.motif = residual.motif;
            summary.signature = residual.signature;
        }
    }

    g_empty_s_motif_summary_cache.emplace(summary_key, summary);
    return summary;
}

bool tryExactEmptySMotifReduction(const VectorRangeTreeMap& start,
                                  const VectorRangeTreeMap& target,
                                  int k,
                                  bool& handled,
                                  bool& result_out) {
    handled = false;
    result_out = false;
    const int motif_m = emptySExactMotifM();
    if (motif_m <= 0) {
        return false;
    }
    if (TreesEqual(start, target)) {
        handled = true;
        result_out = true;
        return true;
    }

    bool solved_by_free_edge = false;
    auto layouts = collectEmptySMotifPeelLayouts(start, target, motif_m, solved_by_free_edge);
    for (const auto& layout : layouts) {
        if (layout.solves_entire_state) {
            if (g_profile.enabled) {
                g_profile.s_empty_motif_free_peels++;
            }
            handled = true;
            result_out = (layout.base_cost <= k);
            return true;
        }
        if (layout.base_cost > k || layout.residual_edges >= countInternalEdges(start)) {
            continue;
        }

        int exact_side_cost = exactEmptySMotifDistance(layout.peeled_start, layout.peeled_target, true);
        if (layout.base_cost + exact_side_cost > k) {
            continue;
        }

        if (g_profile.enabled) {
            if (layout.kind == "common") {
                g_profile.s_empty_motif_common_peels++;
            } else {
                g_profile.s_empty_motif_free_peels++;
            }
        }

        int residual_k = k - layout.base_cost - exact_side_cost;
        if (TreesEqual(layout.residual_start, layout.residual_target)) {
            handled = true;
            result_out = true;
            return true;
        }

        bool residual_handled = false;
        bool residual_result = false;
        tryExactEmptySMotifReduction(layout.residual_start,
                                     layout.residual_target,
                                     residual_k,
                                     residual_handled,
                                     residual_result);
        if (residual_handled && residual_result) {
            handled = true;
            result_out = true;
            return true;
        }
    }

    const int current_edges = countInternalEdges(start);
    const std::string current_motif = classifyEmptySCoreMotif(start, target);
    if (current_edges <= motif_m && motifExactFallbackEligible(current_motif)) {
        int exact_dist = exactEmptySMotifDistance(start, target, true);
        handled = true;
        result_out = (exact_dist <= k);
        return true;
    }

    return false;
}

struct InflightBudgetGuard {
    std::unordered_map<Key128, int, Key128Hash>* map = nullptr;
    Key128 key{};
    bool had_prev = false;
    int prev_budget = -1;
    ~InflightBudgetGuard() {
        if (!map) return;
        if (had_prev) {
            (*map)[key] = prev_budget;
        } else {
            map->erase(key);
        }
    }
};

} // namespace

// Definition: Decide if rotation distance <= k via Li–Xia branching
// Parameters: T_init: source tree; T_final: target tree; k: budget
// Returns: true if a sequence of <= k rotations exists, else false
// Errors: returns false on invalid inputs or internal copy/rotation failures
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k) {
    initProfile();
    if (g_profile.enabled) {
        g_profile.calls_flipdist++;
    }
    ScopedTimer _t(g_profile.enabled ? &g_profile.time_flipdist_ms : nullptr);
    throwIfProfileAbort();
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    const Key128 memo_key = makeKeyFlip(T_init, T_final, k);
    const Key128 pair_key = makeKeyPair(T_init, T_final);
    const Key128 reverse_memo_key = makeKeyFlip(T_final, T_init, k);
    const Key128 reverse_pair_key = makeKeyPair(T_final, T_init);
    auto memo_it = g_memo_flipdist.find(memo_key);
    if (memo_it != g_memo_flipdist.end()) {
        return memo_it->second;
    }
    bool pair_bounds_value = false;
    if (tryBoundsPrune(g_kbounds, pair_key, k, pair_bounds_value)) {
        g_memo_flipdist[memo_key] = pair_bounds_value;
        return pair_bounds_value;
    }
    auto memoReturn = [&](bool value) {
        g_memo_flipdist[memo_key] = value;
        updateBounds(g_kbounds, pair_key, k, value);
        g_memo_flipdist[reverse_memo_key] = value;
        updateBounds(g_kbounds, reverse_pair_key, k, value);
        if (value) {
            updatePairIncumbentUpper(pair_key, k);
            updatePairIncumbentUpper(reverse_pair_key, k);
        }
        return value;
    };

    // BASE CASE: Check if trees are already identical
    if (TreesEqual(T_init, T_final)) {
        debugPrint("Trees already equal, returning true");
        return memoReturn(true);
    }

    static const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> empty_s;
    static const std::vector<std::pair<int, int>> empty_i;
    int conflict_lb = cachedConflictCount(T_init, T_final);
    int pair_lb = heuristicLowerBound(T_init, T_final, HeuristicStateKind::Pair, empty_s, empty_i, conflict_lb);
    if (pair_lb > k) {
        return memoReturn(false);
    }

    bool handled_common = false;
    const bool use_common_edge = tdsCommonEdgeRotatePairEnabled();
    if (use_common_edge) {
        const int min_common_edges = tdsCommonEdgeMinCount();
        const int common_edge_budget_cap = std::max(4, min_common_edges * 5);
        if (k <= common_edge_budget_cap &&
            pairHasCommonEdgesCached(T_init, T_final, min_common_edges)) {
            if (tryCommonEdgeDecomposeCached(T_init, T_final, k, handled_common)) {
                debugPrint("Common edge decomposition succeeded");
                return memoReturn(true);
            }
        }
    }
    if (handled_common) {
        debugPrint("Common edge decomposition failed");
        return memoReturn(false);
    }

    auto [hasFree, freeEdge] = cachedFindFreeEdge(T_init, T_final);
    if (hasFree) {
        debugPrint("FlipDistTree: free edge found, delegating to TreeDistS");
        return memoReturn(TreeDistS(T_init, T_final, k, {}));
    }

    debugPrint("FlipDistTree: budget k=" + std::to_string(k));

    auto edges = getInternalEdges(T_init);
    debugPrint("Found " + std::to_string(edges.size()) + " internal edges");

    std::vector<std::vector<std::pair<int,int>>> independentSubsets;
    std::vector<std::pair<int,int>> current;

    //                       branch on two choices: (1) include e in I; (2) exclude e from I"
    generateAllIndependentSubsets(edges, 0, current, independentSubsets);
    if (g_profile.enabled) {
        g_profile.independent_subsets_initial += (long long)independentSubsets.size();
    }

    debugPrint("Generated " + std::to_string(independentSubsets.size()) + " independent subsets (all)");

    //                       If FlipDistTree-I(T_init, T_final, |I|) returns True, then True"
    for (size_t i = 0; i < independentSubsets.size(); i++) {
        const auto& subset = independentSubsets[i];
        if (subset.empty()) continue;  // Skip empty subsets (I ≠ ∅ requirement)
        if (g_profile.enabled) {
            if ((long long)subset.size() > g_profile.max_i_size) {
                g_profile.max_i_size = (long long)subset.size();
            }
        }

        debugPrint("Trying subset " + std::to_string(i) + " with " + std::to_string(subset.size()) + " edges");

        try {
            // Call TreeDistI (which is FlipDistTree-I from pseudocode)
            if (TreeDistI(T_init, T_final, k, subset)) {
                debugPrint("Found solution with subset " + std::to_string(i));
                return memoReturn(true);
            }
        } catch (const SearchAbortException&) {
            throw;
        } catch (...) {
            debugPrint("Exception in TreeDistI for subset " + std::to_string(i));
            continue;
        }
    }

    debugPrint("No solution found, returning false");
    return memoReturn(false);
}

// Definition: Apply an independent set of rotations, then branch on S
// Parameters: T_init/T_final: trees; k: total budget; I: independent edges
// Returns: true if solvable within k, else false
// Errors: returns false on invalid edges or copy/rotation failures
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I) {
    initProfile();
    if (g_profile.enabled) {
        g_profile.calls_tdi++;
        if ((long long)I.size() > g_profile.max_i_size) {
            g_profile.max_i_size = (long long)I.size();
        }
    }
    ScopedTimer _t(g_profile.enabled ? &g_profile.time_tdi_ms : nullptr);
    throwIfProfileAbort();
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    const Key128 memo_key = makeKeyI(T_init, T_final, k, I);
    const Key128 pair_key = makeKeyPair(T_init, T_final);
    auto memo_it = g_memo_tdi.find(memo_key);
    if (memo_it != g_memo_tdi.end()) {
        return memo_it->second;
    }
    const Key128 bounds_key = makeKeyIBase(T_init, T_final, I);
    bool bounds_value = false;
    if (tryBoundsPrune(g_tdi_bounds, bounds_key, k, bounds_value)) {
        if (g_profile.enabled) g_profile.tdi_bounds_prunes++;
        g_memo_tdi[memo_key] = bounds_value;
        return bounds_value;
    }

    auto inflight_it = g_tdi_inflight_best_budget.find(bounds_key);
    if (tdiInflightPruningEnabled()) {
        if (inflight_it != g_tdi_inflight_best_budget.end() && inflight_it->second >= k) {
            if (g_profile.enabled) g_profile.dominance_prunes++;
            g_memo_tdi[memo_key] = false;
            return false;
        }
        InflightBudgetGuard inflight_guard;
        inflight_guard.map = &g_tdi_inflight_best_budget;
        inflight_guard.key = bounds_key;
        if (inflight_it != g_tdi_inflight_best_budget.end()) {
            inflight_guard.had_prev = true;
            inflight_guard.prev_budget = inflight_it->second;
            if (k > inflight_it->second) {
                g_tdi_inflight_best_budget[bounds_key] = k;
            }
        } else {
            g_tdi_inflight_best_budget[bounds_key] = k;
        }
    } else {
        InflightBudgetGuard inflight_guard;
        (void)inflight_guard;
    }

    auto memoReturn = [&](bool value) {
        g_memo_tdi[memo_key] = value;
        updateBounds(g_tdi_bounds, bounds_key, k, value);
        if (value) {
            updatePairIncumbentUpper(pair_key, k);
        }
        return value;
    };

    // Strong lower bound gate: conflict count + learned failure bounds
    int conflict_lb = cachedConflictCount(T_init, T_final);
    int lower_bound = combinedLowerBoundI(T_init, T_final, I, conflict_lb);
    if (lower_bound > k) {
        if (g_profile.enabled) g_profile.tdi_conflict_lb_prunes++;
        return memoReturn(false);
    }

    int remaining_budget = k - (int)I.size();  // k - |I| from pseudocode
    const bool use_prefix_lb_prune = tdiPrefixPruneEnabled();

    if (remaining_budget < 0) {  // This covers φ(T_init) > k − |I| case
        debugPrint("TreeDistI: Not enough budget for |I| rotations");
        if (g_profile.enabled) g_profile.tdi_remaining_budget_prunes++;
        return memoReturn(false);
    }

    // Special handling when budget exactly equals |I| - try direct solution
    if (remaining_budget == 0) {
        VectorRangeTreeMap T_test = safeCopyTree(T_init);

        // Apply all rotations in I
        for (const auto& edge : I) {
            auto oriented = orientEdge(T_test, edge);
            int parent = oriented.first;
            int child = oriented.second;
            if (parent < 0 || child < 0) continue;

            try {
                if (T_test.getLeftChild(parent) == child) {
                    T_test.rotateRight(parent);
                } else if (T_test.getRightChild(parent) == child) {
                    T_test.rotateLeft(parent);
                }
            } catch (const SearchAbortException&) {
                throw;
            } catch (...) {
                continue;
            }
        }

        if (TreesEqual(T_test, T_final)) {
            debugPrint("TreeDistI: Solved with exactly |I| rotations");
            if (g_profile.enabled) g_profile.tdi_exact_budget_success++;
            return memoReturn(true);
        } else {
            debugPrint("TreeDistI: No budget left and not solved");
            if (g_profile.enabled) g_profile.tdi_exact_budget_fail++;
            return memoReturn(false);
        }
    }

    debugPrint("TreeDistI: Proceeding with remaining_budget=" + std::to_string(remaining_budget));

    if (TreesEqual(T_init, T_final)) {
        debugPrint("TreeDistI: Trees already equal");
        if (g_profile.enabled) g_profile.tdi_tree_equal_success++;
        return memoReturn(true);
    }

    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S;
    std::vector<std::pair<int,int>> newDiagonals;
    std::unordered_set<Key128, Key128Hash> prefix_seen_pairs;
    if (use_prefix_lb_prune) {
        prefix_seen_pairs.reserve(std::max<size_t>(8, I.size() * 2));
    }
    int applied_rotations = 0;

    VectorRangeTreeMap T_bar = safeCopyTree(T_init);  // T̄_init from pseudocode
    if (T_bar.original_nodes.empty()) {
        debugPrint("TreeDistI: Failed to copy tree");
        return memoReturn(false);
    }

    for (const auto& edge : I) {
        auto oriented = orientEdge(T_bar, edge);
        int parent = oriented.first;
        int child = oriented.second;

        debugPrint("TreeDistI: Processing edge (" + std::to_string(parent) + "," + std::to_string(child) + ")");

        if (parent < 0 || child < 0) {
            debugPrint("TreeDistI: Invalid edge, skipping");
            if (g_profile.enabled) g_profile.tdi_rotation_skips++;
            continue;
        }

        int u, v;  // The nodes u,v from pseudocode step 2.2

        if (T_bar.getLeftChild(parent) == child) {
            u = child;   // After rotation, child becomes parent
            v = parent;  // After rotation, parent becomes child
            T_bar.rotateRight(parent);
        } else if (T_bar.getRightChild(parent) == child) {
            u = child;   // After rotation, child becomes parent
            v = parent;  // After rotation, parent becomes child
            T_bar.rotateLeft(parent);
        } else {
            if (g_profile.enabled) g_profile.tdi_rotation_skips++;
            continue;
        }

        debugPrint("TreeDistI: Applied rotation, new edge ē connects " + std::to_string(u) + " and " + std::to_string(v));
        auto cr = T_bar.getRange(v);
        newDiagonals.emplace_back(cr.first, cr.second);
        applied_rotations++;

        if (use_prefix_lb_prune && applied_rotations < (int)I.size()) {
            if (g_profile.enabled) {
                g_profile.tdi_prefix_pairgen_calls++;
            }
            auto pairgen_t0 = g_profile.enabled ? std::chrono::steady_clock::now()
                                                : std::chrono::steady_clock::time_point{};
            appendUniquePartnerPairsFromDiagonal(
                T_bar,
                newDiagonals.back(),
                S,
                prefix_seen_pairs);
            if (g_profile.enabled) {
                g_profile.time_tdi_prefix_pairgen_ms +=
                    std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - pairgen_t0).count();
            }
            const int prefix_remaining_budget = k - applied_rotations;
            if (g_profile.enabled) {
                g_profile.tdi_prefix_conflict_calls++;
            }
            auto conflict_t0 = g_profile.enabled ? std::chrono::steady_clock::now()
                                                 : std::chrono::steady_clock::time_point{};
            const int prefix_conflict_lb = cachedConflictCount(T_bar, T_final);
            if (g_profile.enabled) {
                g_profile.time_tdi_prefix_conflict_ms +=
                    std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - conflict_t0).count();
            }
            if (g_profile.enabled) {
                g_profile.tdi_prefix_lb_checks++;
            }
            if (prefix_conflict_lb > prefix_remaining_budget) {
                if (g_profile.enabled) {
                    g_profile.tdi_prefix_lb_prunes++;
                }
                return memoReturn(false);
            }
            if ((int)S.size() >= tdiPrefixPruneMinSSize()) {
                if (g_profile.enabled) {
                    g_profile.tdi_prefix_bound_calls++;
                }
                auto bound_t0 = g_profile.enabled ? std::chrono::steady_clock::now()
                                                  : std::chrono::steady_clock::time_point{};
                const int prefix_lower_bound = combinedLowerBoundS(T_bar, T_final, S, prefix_conflict_lb);
                if (g_profile.enabled) {
                    g_profile.time_tdi_prefix_bound_ms +=
                        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - bound_t0).count();
                }
                if (prefix_lower_bound > prefix_remaining_budget) {
                    if (g_profile.enabled) {
                        g_profile.tdi_prefix_lb_prunes++;
                    }
                    return memoReturn(false);
                }
            }
        }
    }

    if (!use_prefix_lb_prune) {
        appendPartnerPairsFromDiagonals(T_bar, newDiagonals, S);
    }

    int post_i_conflict_lb = cachedConflictCount(T_bar, T_final);
    int post_i_lower_bound = combinedLowerBoundS(T_bar, T_final, S, post_i_conflict_lb);
    if (post_i_lower_bound > k - (int)I.size()) {
        if (tdiPostIProfileEnabled()) {
            TreeStructureProfile start_profile = describeTreeStructure(T_bar);
            TreeStructureProfile target_profile = describeTreeStructure(T_final);
            noteTdiPostIPruneSample(
                makeKeySBase(T_bar, T_final, S),
                k - (int)I.size(),
                post_i_lower_bound,
                post_i_conflict_lb,
                static_cast<int>(S.size()),
                start_profile.branching_nodes,
                target_profile.branching_nodes,
                start_profile.degree_histogram,
                target_profile.degree_histogram,
                currentTreeToCanonicalString(T_bar),
                currentTreeToCanonicalString(T_final));
        }
        if (g_profile.enabled) g_profile.tdi_post_i_lb_prunes++;
        return memoReturn(false);
    }
    if (g_profile.enabled) {
        if (S.empty()) {
            g_profile.tdi_post_i_empty_s++;
        } else {
            g_profile.tdi_post_i_nonempty_s++;
        }
    }
    bool tds_ok = TreeDistS(T_bar, T_final, k - (int)I.size(), S);
    if (g_profile.enabled) {
        if (tds_ok) g_profile.tdi_tds_success++;
        else g_profile.tdi_tds_fail++;
    }
    return memoReturn(tds_ok);
}

// Definition: Resolve free-edge/partition cases, else branch on S pairs
// Parameters: T_init/T_end: trees; k: budget; S: edge-pair set
// Returns: true if solvable within k, else false
// Errors: returns false on invalid inputs or internal failures
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    initProfile();
    if (g_profile.enabled) {
        g_profile.calls_tds++;
        if ((long long)S.size() > g_profile.max_s_size) {
            g_profile.max_s_size = (long long)S.size();
        }
    }
    ScopedTimer _t(g_profile.enabled ? &g_profile.time_tds_ms : nullptr);
    throwIfProfileAbort();
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_norm;
    if (!S.empty()) {
        S_norm = S;
        dedupeSPairs(S_norm);
    }
    debugPrint("Entering TreeDistS with k=" + std::to_string(k) + ", |S|=" + std::to_string(S_norm.size()));

    const bool s_norm_empty = S_norm.empty();
    const EmptySKeys entry_empty_keys = s_norm_empty ? makeEmptySKeys(T_init, T_end, k) : EmptySKeys{};
    const Key128 memo_key = s_norm_empty ? entry_empty_keys.exact_key : makeKeyS(T_init, T_end, k, S_norm);
    const Key128 pair_key = s_norm_empty ? entry_empty_keys.pair_key : makeKeyPair(T_init, T_end);
    const Key128 bounds_key = s_norm_empty ? entry_empty_keys.bounds_key : makeKeySBase(T_init, T_end, S_norm);
    EmptySFreeHintDisableGuard free_hint_disable_guard(
        s_norm_empty &&
        g_empty_s_disable_free_hint_depth == 0 &&
        shouldDisableEmptySFreeHintsForRootShape(T_init, T_end, pair_key));
    auto memo_it = g_memo_tds.find(memo_key);
    if (memo_it != g_memo_tds.end()) {
        if (g_profile.enabled) {
            g_profile.memo_hits_tds++;
        }
        const bool memo_value = memo_it->second;
        updateBounds(g_tds_bounds, bounds_key, k, memo_value);
        if (memo_value) {
            updatePairIncumbentUpper(pair_key, k);
        }
        return memo_value;
    }
    bool bounds_value = false;
    if (tryBoundsPrune(g_tds_bounds, bounds_key, k, bounds_value)) {
        if (g_profile.enabled) {
            g_profile.bounds_hits_tds++;
        }
        g_memo_tds[memo_key] = bounds_value;
        if (bounds_value) {
            updatePairIncumbentUpper(pair_key, k);
        }
        return bounds_value;
    }
    bool pair_known_value = false;
    if (tryBoundsPrune(g_kbounds, pair_key, k, pair_known_value)) {
        if (!pair_known_value) {
            g_memo_tds[memo_key] = false;
            updateBounds(g_tds_bounds, bounds_key, k, false);
            return false;
        }
        if (s_norm_empty) {
            g_memo_tds[memo_key] = true;
            updateBounds(g_tds_bounds, bounds_key, k, true);
            updatePairIncumbentUpper(pair_key, k);
            return true;
        }
    }

    auto inflight_it = g_tds_inflight_best_budget.find(bounds_key);
    if (tdsInflightPruningEnabled() && inflight_it != g_tds_inflight_best_budget.end() && inflight_it->second >= k) {
        if (g_profile.enabled) g_profile.dominance_prunes++;
        g_memo_tds[memo_key] = false;
        return false;
    }
    InflightBudgetGuard inflight_guard;
    if (tdsInflightPruningEnabled()) {
        inflight_guard.map = &g_tds_inflight_best_budget;
        inflight_guard.key = bounds_key;
        if (inflight_it != g_tds_inflight_best_budget.end()) {
            inflight_guard.had_prev = true;
            inflight_guard.prev_budget = inflight_it->second;
            if (k > inflight_it->second) {
                g_tds_inflight_best_budget[bounds_key] = k;
            }
        } else {
            g_tds_inflight_best_budget[bounds_key] = k;
        }
    }

    auto memoReturn = [&](bool value) {
        g_memo_tds[memo_key] = value;
        updateBounds(g_tds_bounds, bounds_key, k, value);
        if (value) {
            updatePairIncumbentUpper(pair_key, k);
        }
        return value;
    };

    if (k < 0) {
        debugPrint("TreeDistS: Negative budget");
        return memoReturn(false);
    }
    int conflict_lb = cachedConflictCountWithKey(T_init, T_end, pair_key);
    if (conflict_lb > k) {
        return memoReturn(false);
    }
    if (conflict_lb == 0 && TreesEqual(T_init, T_end)) {
        debugPrint("TreeDistS: Trees already equal");
        return memoReturn(true);
    }
    if (k == 0) {
        debugPrint("TreeDistS: No budget and trees differ");
        return memoReturn(false);
    }

    if (s_norm_empty && pairIncumbentUpper(pair_key) <= k) {
        return memoReturn(true);
    }

    int lower_bound = (s_norm_empty && defaultEmptySBoundConfig())
                          ? emptySLowerBoundFromKeys(conflict_lb, pair_key, bounds_key)
                          : combinedLowerBoundS(T_init, T_end, S_norm, conflict_lb);
    if (lower_bound > k) {
        return memoReturn(false);
    }
    const int incumbent_upper = pairIncumbentUpper(pair_key);
    if (k < incumbent_upper && lower_bound >= incumbent_upper) {
        return memoReturn(false);
    }

    VectorRangeTreeMap T_work = safeCopyTree(T_init);
    int k_work = k;
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_work = S_norm;
    struct UnfoldedStateMemoEntry {
        Key128 exact_key;
        Key128 bounds_key;
        int budget = 0;
    };
    std::vector<UnfoldedStateMemoEntry> unfolded_bounds_entries;

    auto commitResult = [&](bool value) {
        for (const auto& entry : unfolded_bounds_entries) {
            g_memo_tds[entry.exact_key] = value;
            updateBounds(g_tds_bounds, entry.bounds_key, entry.budget, value);
        }
        return memoReturn(value);
    };

    Key128 work_pair_key{};
    while (k_work > 0) {
        EmptySKeys work_empty_keys{};
        const bool s_work_empty = S_work.empty();
        if (s_work_empty) {
            work_empty_keys = makeEmptySKeys(T_work, T_end, k_work);
            work_pair_key = work_empty_keys.pair_key;
        } else {
            work_pair_key = makeKeyPair(T_work, T_end);
        }
        if (s_work_empty && pairIncumbentUpper(work_pair_key) <= k_work) {
            return commitResult(true);
        }

        int work_conflict_lb = cachedConflictCountWithKey(T_work, T_end, work_pair_key);
        int work_lb = (s_work_empty && defaultEmptySBoundConfig())
                          ? emptySLowerBoundFromKeys(work_conflict_lb,
                                                     work_pair_key,
                                                     work_empty_keys.bounds_key)
                          : combinedLowerBoundS(T_work, T_end, S_work, work_conflict_lb);
        if (work_lb > k_work) {
            return commitResult(false);
        }

        const Key128 state_key = s_work_empty ? work_empty_keys.exact_key : makeKeyS(T_work, T_end, k_work, S_work);
        auto state_it = g_memo_tds.find(state_key);
        if (state_it != g_memo_tds.end()) {
            if (g_profile.enabled) {
                g_profile.memo_hits_tds++;
            }
            return commitResult(state_it->second);
        }
        const Key128 state_bounds_key = s_work_empty ? work_empty_keys.bounds_key : makeKeySBase(T_work, T_end, S_work);
        bool state_bounds_value = false;
        if (tryBoundsPrune(g_tds_bounds, state_bounds_key, k_work, state_bounds_value)) {
            if (g_profile.enabled) {
                g_profile.bounds_hits_tds++;
            }
            return commitResult(state_bounds_value);
        }
        unfolded_bounds_entries.push_back({state_key, state_bounds_key, k_work});

        const bool partition_cache_enabled = tdsPartitionCacheEnabled();

        auto [hasFree, freeEdge] = cachedFindFreeEdge(T_work, T_end, work_pair_key);
        if (g_profile.enabled) {
            if (hasFree) g_profile.free_edge_hits++;
            else g_profile.free_edge_misses++;
        }
        if (!hasFree) {
            break;
        }

        debugPrint("TreeDistS: Found free edge (" + std::to_string(freeEdge.first) + "," + std::to_string(freeEdge.second) + ")");

        try {
            int parent = freeEdge.first;
            int child = freeEdge.second;

            std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_filtered;
            {
                ScopedTimer _tfilter(g_profile.enabled ? &g_profile.time_tds_free_filter_ms : nullptr);
                S_filtered.reserve(S_work.size() + 2);
                for (const auto& pair : S_work) {
                    if (!(pair.first == freeEdge || pair.second == freeEdge)) {
                        S_filtered.push_back(pair);
                    }
                }
            }
            debugPrint("TreeDistS: Filtered S from " + std::to_string(S_work.size()) + " to " + std::to_string(S_filtered.size()) + " pairs");

            if (!refreshOriginalPreorderInPlace(T_work)) {
                return commitResult(false);
            }
            VectorRangeTreeMap& T_bar = T_work;
            int u, v;

            {
                ScopedTimer _trotate(g_profile.enabled ? &g_profile.time_tds_free_rotate_ms : nullptr);
                if (T_bar.getLeftChild(parent) == child) {
                    T_bar.rotateRight(parent);
                    u = child;
                    v = parent;
                } else if (T_bar.getRightChild(parent) == child) {
                    T_bar.rotateLeft(parent);
                    u = child;
                    v = parent;
                } else {
                    return commitResult(false);
                }
            }

            if (TreesEqual(T_bar, T_end)) {
                debugPrint("TreeDistS: Solved with free edge rotation");
                return commitResult(true);
            }

            auto cr = T_bar.getRange(v);
            int L = cr.first;
            int R = cr.second;

            {
                ScopedTimer _tpartner(g_profile.enabled ? &g_profile.time_tds_free_partner_ms : nullptr);
                appendPartnerPairsFromSingleDiagonalProfiled(T_bar, {L, R}, S_filtered, nullptr);
                dedupeSPairs(S_filtered);
            }

            const auto partition_parent_range = T_bar.getRange(u);
            const auto partition_child_range = T_bar.getRange(v);
            const bool precheck_enabled = tdsPartitionPrecheckEnabled();
            const bool split_cache_enabled = tdsPartitionSplitCacheEnabled();
            Key128 partition_cache_pair_key{};
            FreeEdgePartitionStructureKey structure_cache_key{};
            const bool partition_identity_key_needed = partition_cache_enabled || split_cache_enabled;
            if (partition_identity_key_needed) {
                partition_cache_pair_key = makeKeyPair(T_bar, T_end);
                structure_cache_key = FreeEdgePartitionStructureKey{
                    partition_cache_pair_key,
                    static_cast<std::uint32_t>(partition_parent_range.first),
                    static_cast<std::uint32_t>(partition_parent_range.second),
                    static_cast<std::uint32_t>(partition_child_range.first),
                    static_cast<std::uint32_t>(partition_child_range.second)};
            }
            const int partition_budget = k_work - 1;

            if (g_profile.enabled) {
                g_profile.calls_partition++;
            }
            FreeEdgePartitionStructureCacheEntry partition_structure_cache_local;
            FreeEdgePartitionStructureCacheEntry* partition_structure_cache_ptr = &partition_structure_cache_local;
            if (partition_cache_enabled) {
                partition_structure_cache_ptr = &g_free_edge_partition_structure_cache
                    .try_emplace(structure_cache_key, FreeEdgePartitionStructureCacheEntry{})
                    .first->second;
            }
            auto& partition_structure_cache = *partition_structure_cache_ptr;

            if (!partition_structure_cache.computed) {
                if (precheck_enabled) {
                    if (g_profile.enabled) {
                        g_profile.partition_precheck_calls++;
                    }
                }
                ScopedTimer _tp(g_profile.enabled ? &g_profile.time_partition_ms : nullptr);
                try {
                    if (precheck_enabled) {
                        ScopedTimer _tpre(g_profile.enabled ? &g_profile.time_partition_precheck_ms : nullptr);
                        partition_structure_cache.precheck_ok =
                            partitionNodeSetsMatchByChildRange(T_bar, T_end, partition_child_range);
                    } else {
                        partition_structure_cache.precheck_ok = true;
                    }
                    if (g_profile.enabled) {
                        if (partition_structure_cache.precheck_ok) {
                            g_profile.partition_precheck_matches++;
                        } else {
                            g_profile.partition_precheck_mismatches++;
                        }
                    }
                    if (!partition_structure_cache.precheck_ok) {
                        partition_structure_cache.computed = true;
                    } else {
                        ScopedTimer _tbuild(g_profile.enabled ? &g_profile.time_partition_build_ms : nullptr);
                            std::tie(partition_structure_cache.T_bar1, partition_structure_cache.T_bar2) =
                                VectorRangeTreeMap::partitionAlongSubtreeRange(T_bar,
                                                                               partition_child_range);
                        const TargetPartitionRangeCacheKey target_partition_key{
                            makeKeyPair(T_end, T_end),
                            static_cast<std::uint32_t>(partition_child_range.first),
                            static_cast<std::uint32_t>(partition_child_range.second)};
                        auto& target_partition_entry = g_target_partition_range_cache
                            .try_emplace(target_partition_key, TargetPartitionRangeCacheEntry{})
                            .first->second;
                        if (!target_partition_entry.computed) {
                            std::tie(target_partition_entry.T_end1, target_partition_entry.T_end2) =
                                VectorRangeTreeMap::partitionAlongEdge(T_end,
                                                                       partition_parent_range,
                                                                       partition_child_range);
                            target_partition_entry.computed = true;
                        }
                        partition_structure_cache.T_end1_ptr = &target_partition_entry.T_end1;
                        partition_structure_cache.T_end2_ptr = &target_partition_entry.T_end2;

                        if (partition_structure_cache.T_bar1.original_nodes != target_partition_entry.T_end1.original_nodes ||
                            partition_structure_cache.T_bar2.original_nodes != target_partition_entry.T_end2.original_nodes) {
                            partition_structure_cache.partition_nodes_match = false;
                        } else {
                            ScopedTimer _tedges(g_profile.enabled ? &g_profile.time_partition_count_edges_ms : nullptr);
                            partition_structure_cache.n1 = countInternalEdges(partition_structure_cache.T_bar1);
                            partition_structure_cache.n2 = countInternalEdges(partition_structure_cache.T_bar2);
                        }
                        partition_structure_cache.computed = true;
                    }
                } catch (const SearchAbortException&) {
                    throw;
                } catch (...) {
                    partition_structure_cache.had_exception = true;
                    partition_structure_cache.computed = true;
                }
            } else if (precheck_enabled && g_profile.enabled) {
                g_profile.partition_precheck_calls++;
                if (partition_structure_cache.precheck_ok) {
                    g_profile.partition_precheck_matches++;
                } else {
                    g_profile.partition_precheck_mismatches++;
                }
            }

            bool side_cache_key_ready = false;
            FreeEdgePartitionSideCacheKey side_cache_key{};
            FreeEdgePartitionSplitDecisionCacheKey split_cache_key{};
            FreeEdgePartitionSplitSignatureCacheKey split_signature_cache_key{};
            FreeEdgePartitionSideCacheEntry partition_cache_local;
            FreeEdgePartitionSideCacheEntry* partition_cache_ptr = &partition_cache_local;
            std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> side_pairs_left;
            std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> side_pairs_right;

            if (!partition_structure_cache.had_exception &&
                partition_structure_cache.precheck_ok &&
                partition_structure_cache.partition_nodes_match) {
                try {
                    if (partition_cache_enabled || split_cache_enabled) {
                        Key128 partition_side_input_key = makeKeySBase(T_bar, T_end, S_filtered);
                        const Key128 partition_side_signature_key =
                            makePartitionFilteredSignature(partition_cache_pair_key, partition_side_input_key);
                        side_cache_key = FreeEdgePartitionSideCacheKey{
                            structure_cache_key,
                            partition_side_signature_key};
                    }
                    side_cache_key_ready = true;
                    split_cache_key = FreeEdgePartitionSplitDecisionCacheKey{
                        side_cache_key,
                        partition_budget};

                    if (split_cache_enabled) {
                        auto& split_state_cache = g_free_edge_partition_split_budget_cache
                            .try_emplace(side_cache_key, FreeEdgePartitionSplitBudgetCacheEntry{})
                            .first->second;
                        if (split_state_cache.computed) {
                            if (split_state_cache.solved_from >= 0 && partition_budget >= split_state_cache.solved_from) {
                                if (g_profile.enabled) {
                                    g_profile.tds_partition_split_cache_hits++;
                                }
                                return commitResult(true);
                            }
                            if (split_state_cache.failed_up_to >= 0 && partition_budget <= split_state_cache.failed_up_to) {
                                if (g_profile.enabled) {
                                    g_profile.tds_partition_split_cache_hits++;
                                }
                                return commitResult(false);
                            }
                        }
                        const auto split_cache_it_fast = g_free_edge_partition_split_cache.find(split_cache_key);
                        if (split_cache_it_fast != g_free_edge_partition_split_cache.end() &&
                            split_cache_it_fast->second.computed) {
                            if (g_profile.enabled) {
                                g_profile.tds_partition_split_cache_hits++;
                            }
                            split_state_cache.computed = true;
                            if (split_cache_it_fast->second.solvable) {
                                if (split_state_cache.solved_from < 0 || partition_budget < split_state_cache.solved_from) {
                                    split_state_cache.solved_from = partition_budget;
                                    split_state_cache.witness_k1 = split_cache_it_fast->second.witness_k1;
                                }
                                if (split_state_cache.failed_up_to >= split_state_cache.solved_from) {
                                    split_state_cache.failed_up_to = split_state_cache.solved_from - 1;
                                }
                            } else {
                                split_state_cache.failed_up_to = std::max(split_state_cache.failed_up_to, partition_budget);
                            }
                            return commitResult(split_cache_it_fast->second.solvable);
                        }
                        if (g_profile.enabled) {
                            g_profile.tds_partition_split_cache_misses++;
                        }
                    }

                    if (partition_cache_enabled) {
                        partition_cache_ptr = &g_free_edge_partition_side_cache
                            .try_emplace(side_cache_key, FreeEdgePartitionSideCacheEntry{})
                            .first->second;
                    }

                    {
                        ScopedTimer _tsplit(g_profile.enabled ? &g_profile.time_partition_split_s_ms : nullptr);
                        std::tie(side_pairs_left, side_pairs_right) =
                            partitionS(S_filtered, partition_structure_cache.T_bar1, partition_structure_cache.T_bar2);
                    }
                } catch (const SearchAbortException&) {
                    throw;
                } catch (...) {
                        partition_structure_cache.had_exception = true;
                }
            }
            auto& partition_cache = *partition_cache_ptr;

            if (!partition_cache.computed && side_cache_key_ready &&
                (!partition_cache_enabled || partition_cache_ptr != &partition_cache_local)) {
                try {
                    {
                        ScopedTimer _tstats(g_profile.enabled ? &g_profile.time_partition_side_stats_ms : nullptr);
                        partition_cache.S1 = std::move(side_pairs_left);
                        partition_cache.S2 = std::move(side_pairs_right);
                    }
                    {
                        ScopedTimer _tconf(g_profile.enabled ? &g_profile.time_partition_conflicts_ms : nullptr);
                        partition_cache.conf1 = cachedConflictCount(partition_structure_cache.T_bar1,
                                                                  partition_structure_cache.T_end1);
                        partition_cache.conf2 = cachedConflictCount(partition_structure_cache.T_bar2,
                                                                  partition_structure_cache.T_end2);
                    }
                    {
                        ScopedTimer _tbounds(g_profile.enabled ? &g_profile.time_partition_bounds_ms : nullptr);
                        partition_cache.s1_pair_key = makeKeyPair(partition_structure_cache.T_bar1,
                                                                 partition_structure_cache.T_end1);
                        partition_cache.s2_pair_key = makeKeyPair(partition_structure_cache.T_bar2,
                                                                 partition_structure_cache.T_end2);
                        partition_cache.s1_bounds_key = makeKeySBase(partition_structure_cache.T_bar1,
                                                                     partition_structure_cache.T_end1,
                                                                     partition_cache.S1);
                        partition_cache.s2_bounds_key = makeKeySBase(partition_structure_cache.T_bar2,
                                                                     partition_structure_cache.T_end2,
                                                                     partition_cache.S2);
                    }
                    partition_cache.computed = true;
                } catch (const SearchAbortException&) {
                    throw;
                } catch (...) {
                    partition_cache.had_exception = true;
                    partition_cache.computed = true;
                }
            }

            const int snapshot_n1 = partition_structure_cache.n1;
            const int snapshot_n2 = partition_structure_cache.n2;
            const bool snapshot_precheck_ok = partition_structure_cache.precheck_ok;
            const bool snapshot_partition_nodes_match = partition_structure_cache.partition_nodes_match;
            const bool snapshot_structure_exception = partition_structure_cache.had_exception;
            const bool snapshot_side_exception = partition_cache.had_exception;
            const int snapshot_conf1 = partition_cache.conf1;
            const int snapshot_conf2 = partition_cache.conf2;
            const Key128 snapshot_s1_bounds_key = partition_cache.s1_bounds_key;
            const Key128 snapshot_s2_bounds_key = partition_cache.s2_bounds_key;
            const Key128 snapshot_s1_pair_key = partition_cache.s1_pair_key;
            const Key128 snapshot_s2_pair_key = partition_cache.s2_pair_key;
            const VectorRangeTreeMap& snapshot_bar1 = partition_structure_cache.T_bar1;
            const VectorRangeTreeMap& snapshot_bar2 = partition_structure_cache.T_bar2;
            const VectorRangeTreeMap& snapshot_end1 = partition_structure_cache.T_end1_ptr
                                                          ? *partition_structure_cache.T_end1_ptr
                                                          : partition_structure_cache.T_end1;
            const VectorRangeTreeMap& snapshot_end2 = partition_structure_cache.T_end2_ptr
                                                          ? *partition_structure_cache.T_end2_ptr
                                                          : partition_structure_cache.T_end2;
            const auto& snapshot_S1 = partition_cache.S1;
            const auto& snapshot_S2 = partition_cache.S2;
            const bool snapshot_s1_empty = snapshot_S1.empty();
            const bool snapshot_s2_empty = snapshot_S2.empty();
            const bool use_incumbent_prune = incumbentPruneEnabled();
            const int snapshot_s1_pair_incumbent = use_incumbent_prune
                                                      ? pairIncumbentUpper(snapshot_s1_pair_key)
                                                      : std::numeric_limits<int>::max();
            const int snapshot_s2_pair_incumbent = use_incumbent_prune
                                                      ? pairIncumbentUpper(snapshot_s2_pair_key)
                                                      : std::numeric_limits<int>::max();
            const auto boundsMaxFail = [](const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                        const Key128& key) -> int {
                auto it = bounds_map.find(key);
                if (it == bounds_map.end()) return -1;
                return it->second.max_fail;
            };
            const auto pairMaxFail = [](const Key128& key) -> int {
                auto it = g_kbounds.find(key);
                if (it == g_kbounds.end()) return -1;
                return it->second.max_fail;
            };
            const auto sideLowerBound = [](const Key128& pair_key, const Key128& bounds_key, int conflict) -> int {
                const int bounds_lb = requiredBudgetFromBounds(g_tds_bounds, bounds_key);
                const int kbounds_lb = requiredBudgetFromBounds(g_kbounds, pair_key);
                return std::max(conflict, std::max(bounds_lb, kbounds_lb));
            };
            const auto sideKnownFailMax = [&](const Key128& pair_key, const Key128& bounds_key) -> int {
                const int tb_fail = boundsMaxFail(g_tds_bounds, bounds_key);
                const int pair_fail = pairMaxFail(pair_key);
                return pair_fail > tb_fail ? pair_fail : tb_fail;
            };
            const auto exactMinSuccess = [](const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                           const Key128& key) -> int {
                auto it = bounds_map.find(key);
                if (it == bounds_map.end()) return -1;
                if (it->second.min_success < 0) return -1;
                if (it->second.max_fail == it->second.min_success - 1) return it->second.min_success;
                return -1;
            };
            const auto exactMinSuccessOrZero = [&](const Key128& pair_key, const Key128& bounds_key) -> int {
                const int state_success = exactMinSuccess(g_tds_bounds, bounds_key);
                const int pair_success = exactMinSuccess(g_kbounds, pair_key);
                return std::max(state_success, pair_success);
            };
            auto makeSideBudgetCacheKey = [](const Key128& side_pair_key,
                                            const Key128& side_bounds_key) {
                return FreeEdgePartitionSideBudgetCacheKey{
                    side_pair_key,
                    side_bounds_key};
            };
            const auto checkSidePartitionFeasibility = [&](const VectorRangeTreeMap& side_bar,
                                                           const VectorRangeTreeMap& side_end,
                                                           const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>&
                                                                   side_pairs,
                                                           const Key128 side_bounds_key,
                                                           const Key128 side_pair_key,
                                                           const int side_pair_incumbent,
                                                           const bool side_empty,
                                                           const int side_budget,
                                                           const bool is_side1) -> bool {
                if (side_empty && side_pair_incumbent < std::numeric_limits<int>::max() &&
                    side_pair_incumbent <= side_budget) {
                    auto& state_cache_entry = g_free_edge_partition_side_budget_cache
                                                 .try_emplace(makeSideBudgetCacheKey(side_pair_key, side_bounds_key),
                                                              FreeEdgePartitionSideBudgetCacheEntry{})
                                                 .first->second;
                    state_cache_entry.computed = true;
                    state_cache_entry.known_min_feasible_budget =
                        std::min(state_cache_entry.known_min_feasible_budget, side_budget);
                    return true;
                }

                bool pair_known = false;
                if (tryBoundsPrune(g_kbounds, side_pair_key, side_budget, pair_known)) {
                    if (!pair_known) {
                        if (g_profile.enabled) {
                            if (is_side1) {
                                g_profile.partition_side1_bound_hits++;
                            } else {
                                g_profile.partition_side2_bound_hits++;
                            }
                        }
                        auto& side_cache_entry = g_free_edge_partition_side_budget_cache
                                                     .try_emplace(makeSideBudgetCacheKey(side_pair_key, side_bounds_key),
                                                                  FreeEdgePartitionSideBudgetCacheEntry{})
                                                     .first->second;
                        side_cache_entry.computed = true;
                        side_cache_entry.known_not_feasible_below =
                            std::max(side_cache_entry.known_not_feasible_below, side_budget);
                        return false;
                    }
                    if (side_empty) {
                        auto& side_cache_entry = g_free_edge_partition_side_budget_cache
                                                     .try_emplace(makeSideBudgetCacheKey(side_pair_key, side_bounds_key),
                                                                  FreeEdgePartitionSideBudgetCacheEntry{})
                                                     .first->second;
                        side_cache_entry.computed = true;
                        side_cache_entry.known_min_feasible_budget =
                            std::min(side_cache_entry.known_min_feasible_budget, side_budget);
                        return true;
                    }
                }

                bool known = false;
                const bool known_prune = tryBoundsPrune(g_tds_bounds, side_bounds_key, side_budget, known);
                if (known_prune) {
                    if (g_profile.enabled) {
                        if (is_side1) {
                            g_profile.partition_side1_bound_hits++;
                        } else {
                            g_profile.partition_side2_bound_hits++;
                        }
                    }
                    auto& side_cache_entry = g_free_edge_partition_side_budget_cache
                                                 .try_emplace(makeSideBudgetCacheKey(side_pair_key, side_bounds_key),
                                                              FreeEdgePartitionSideBudgetCacheEntry{})
                                                 .first->second;
                    side_cache_entry.computed = true;
                    if (known) {
                        side_cache_entry.known_min_feasible_budget =
                            std::min(side_cache_entry.known_min_feasible_budget, side_budget);
                    } else {
                        side_cache_entry.known_not_feasible_below =
                            std::max(side_cache_entry.known_not_feasible_below, side_budget);
                    }
                    return known;
                }

                const auto side_cache_key = makeSideBudgetCacheKey(side_pair_key, side_bounds_key);
                const auto side_cache_it = g_free_edge_partition_side_budget_cache.find(side_cache_key);
                if (side_cache_it != g_free_edge_partition_side_budget_cache.end() &&
                    side_cache_it->second.computed) {
                    auto& side_cache_entry = side_cache_it->second;
                    if (side_cache_entry.known_not_feasible_below >= side_budget) {
                        if (g_profile.enabled) {
                            g_profile.tds_partition_side_budget_cache_hits++;
                        }
                        return false;
                    }
                    if (side_cache_entry.known_min_feasible_budget >= 0 &&
                        side_cache_entry.known_min_feasible_budget <= side_budget) {
                        if (g_profile.enabled) {
                            g_profile.tds_partition_side_budget_cache_hits++;
                        }
                        return true;
                    }
                }

                const bool feasible = TreeDistS(side_bar, side_end, side_budget, side_pairs);
                if (g_profile.enabled) {
                    if (is_side1) {
                        g_profile.partition_side1_recursions++;
                    } else {
                        g_profile.partition_side2_recursions++;
                    }
                }
                auto& side_cache_entry = g_free_edge_partition_side_budget_cache
                                             .try_emplace(side_cache_key, FreeEdgePartitionSideBudgetCacheEntry{})
                                             .first->second;
                side_cache_entry.computed = true;
                if (feasible) {
                    side_cache_entry.known_min_feasible_budget =
                        std::min(side_cache_entry.known_min_feasible_budget, side_budget);
                } else {
                    side_cache_entry.known_not_feasible_below =
                        std::max(side_cache_entry.known_not_feasible_below, side_budget);
                }
                return feasible;
            };

            const auto commitSplitDecision = [&](bool solvable, int witness_k1 = -1) -> bool {
                if (split_cache_enabled) {
                    auto& split_entry = g_free_edge_partition_split_cache
                                            .try_emplace(split_cache_key, FreeEdgePartitionSplitDecisionCacheEntry{})
                                            .first->second;
                    split_entry.computed = true;
                    split_entry.solvable = solvable;
                    split_entry.had_exception = false;
                    split_entry.witness_k1 = solvable ? witness_k1 : -1;

                    auto& split_state_cache = g_free_edge_partition_split_budget_cache
                                                 .try_emplace(split_cache_key.side_cache_key,
                                                              FreeEdgePartitionSplitBudgetCacheEntry{})
                                                 .first->second;
                    split_state_cache.computed = true;
                    if (solvable) {
                        if (split_state_cache.solved_from < 0 || partition_budget < split_state_cache.solved_from) {
                            split_state_cache.solved_from = partition_budget;
                            split_state_cache.witness_k1 = witness_k1;
                            if (split_state_cache.failed_up_to >= split_state_cache.solved_from) {
                                split_state_cache.failed_up_to = split_state_cache.solved_from - 1;
                            }
                        }
                    } else {
                        split_state_cache.failed_up_to = std::max(split_state_cache.failed_up_to, partition_budget);
                    }
                    if (split_cache_enabled) {
                        auto& split_signature_cache =
                            g_free_edge_partition_split_signature_cache
                                .try_emplace(split_signature_cache_key,
                                             FreeEdgePartitionSplitBudgetCacheEntry{})
                                .first->second;
                        split_signature_cache.computed = true;
                        if (solvable) {
                            if (split_signature_cache.solved_from < 0 ||
                                partition_budget < split_signature_cache.solved_from) {
                                split_signature_cache.solved_from = partition_budget;
                                split_signature_cache.witness_k1 = witness_k1;
                                if (split_signature_cache.failed_up_to >= split_signature_cache.solved_from) {
                                    split_signature_cache.failed_up_to = split_signature_cache.solved_from - 1;
                                }
                            }
                        } else {
                            split_signature_cache.failed_up_to = std::max(split_signature_cache.failed_up_to,
                                                                          partition_budget);
                        }
                    }
                }
                return commitResult(solvable);
            };

            if (split_cache_enabled && side_cache_key_ready && !snapshot_structure_exception &&
                !snapshot_side_exception && snapshot_precheck_ok && snapshot_partition_nodes_match) {
                split_signature_cache_key = makeSplitSignatureCacheKey(snapshot_s1_bounds_key,
                                                                       snapshot_s2_bounds_key);
                auto split_signature_cache_it = g_free_edge_partition_split_signature_cache.find(split_signature_cache_key);
                    if (split_signature_cache_it != g_free_edge_partition_split_signature_cache.end()) {
                        auto& split_signature_cache = split_signature_cache_it->second;
                        if (split_signature_cache.computed) {
                        if (split_signature_cache.solved_from >= 0 &&
                            partition_budget >= split_signature_cache.solved_from) {
                            if (g_profile.enabled) {
                                g_profile.tds_partition_split_cache_hits++;
                            }
                            return commitSplitDecision(true, split_signature_cache.witness_k1);
                        }
                        if (split_signature_cache.failed_up_to >= 0 &&
                            partition_budget <= split_signature_cache.failed_up_to) {
                            if (g_profile.enabled) {
                                g_profile.tds_partition_split_cache_hits++;
                            }
                            return commitSplitDecision(false);
                        }
                    }
                }
            }

            if (snapshot_structure_exception || snapshot_side_exception) {
                debugPrint("TreeDistS: Partitioning failed, continuing free-edge contraction");
                S_work = std::move(S_filtered);
                --k_work;
                continue;
            }
            if (!snapshot_precheck_ok || !snapshot_partition_nodes_match) {
                debugPrint("TreeDistS: Partition mismatch, continuing free-edge contraction");
                S_work = std::move(S_filtered);
                --k_work;
                continue;
            }

            if (snapshot_n1 == 0) {
                if (g_profile.enabled) g_profile.partition_empty_side_shortcuts++;
                if (snapshot_conf2 > k_work - 1) {
                    return commitSplitDecision(false);
                }
                const bool side2_ok = checkSidePartitionFeasibility(snapshot_bar2,
                                                                   snapshot_end2,
                                                                   snapshot_S2,
                                                                   snapshot_s2_bounds_key,
                                                                   snapshot_s2_pair_key,
                                                                   snapshot_s2_pair_incumbent,
                                                                   snapshot_s2_empty,
                                                                   partition_budget,
                                                                   false);
                return commitSplitDecision(side2_ok, 0);
            }

            if (snapshot_n2 == 0) {
                if (g_profile.enabled) g_profile.partition_empty_side_shortcuts++;
                if (snapshot_conf1 > k_work - 1) {
                    return commitSplitDecision(false);
                }
                const bool side1_ok = checkSidePartitionFeasibility(snapshot_bar1,
                                                                   snapshot_end1,
                                                                   snapshot_S1,
                                                                   snapshot_s1_bounds_key,
                                                                   snapshot_s1_pair_key,
                                                                   snapshot_s1_pair_incumbent,
                                                                   snapshot_s1_empty,
                                                                   partition_budget,
                                                                   true);
                return commitSplitDecision(side1_ok, partition_budget);
            }

            const int side1_fail_max = sideKnownFailMax(snapshot_s1_pair_key, snapshot_s1_bounds_key);
            const int side2_fail_max = sideKnownFailMax(snapshot_s2_pair_key, snapshot_s2_bounds_key);
            const int side1_state_exact_success = exactMinSuccess(g_tds_bounds, snapshot_s1_bounds_key);
            const int side2_state_exact_success = exactMinSuccess(g_tds_bounds, snapshot_s2_bounds_key);
            const int side1_pair_exact_success = exactMinSuccess(g_kbounds, snapshot_s1_pair_key);
            const int side2_pair_exact_success = exactMinSuccess(g_kbounds, snapshot_s2_pair_key);
            const int side1_lb = sideLowerBound(snapshot_s1_pair_key, snapshot_s1_bounds_key, snapshot_conf1);
            const int side2_lb = sideLowerBound(snapshot_s2_pair_key, snapshot_s2_bounds_key, snapshot_conf2);
            const int side1_exact_lb = std::max(side1_state_exact_success, side1_pair_exact_success);
            const int side2_exact_lb = std::max(side2_state_exact_success, side2_pair_exact_success);
            int partition_min_k1 = std::max(0, side1_lb);
            if (side1_fail_max >= 0) {
                partition_min_k1 = std::max(partition_min_k1, side1_fail_max + 1);
            }
            if (side1_exact_lb >= 0) {
                partition_min_k1 = std::max(partition_min_k1, side1_exact_lb);
            }

            int partition_max_k1 = std::min(partition_budget, partition_budget - side2_lb);
            if (side2_fail_max >= 0) {
                partition_max_k1 = std::min(partition_max_k1, partition_budget - (side2_fail_max + 1));
            }
            if (side2_exact_lb >= 0) {
                partition_max_k1 = std::min(partition_max_k1, partition_budget - side2_exact_lb);
            }
            partition_max_k1 = std::min(partition_max_k1, partition_budget);
            if (side1_lb + side2_lb > partition_budget) {
                return commitSplitDecision(false);
            }
            if (side1_exact_lb >= 0 && side2_exact_lb >= 0) {
                const int guaranteed_k1_min = std::max(partition_min_k1, side1_exact_lb);
                const int guaranteed_k1_max = std::min(partition_max_k1, partition_budget - side2_exact_lb);
                if (guaranteed_k1_min <= guaranteed_k1_max) {
                    return commitSplitDecision(true, guaranteed_k1_min);
                }
            }
            if (snapshot_s1_empty && snapshot_s1_pair_incumbent < std::numeric_limits<int>::max() &&
                side2_exact_lb >= 0) {
                const int guaranteed_k1_min = std::max(partition_min_k1, snapshot_s1_pair_incumbent);
                const int guaranteed_k1_max = std::min(partition_max_k1, partition_budget - side2_exact_lb);
                if (guaranteed_k1_min <= guaranteed_k1_max) {
                    return commitSplitDecision(true, guaranteed_k1_min);
                }
            }
            if (snapshot_s2_empty && snapshot_s2_pair_incumbent < std::numeric_limits<int>::max() &&
                side1_exact_lb >= 0) {
                const int guaranteed_k1_min = std::max(partition_min_k1, side1_exact_lb);
                const int guaranteed_k1_max = std::min(partition_max_k1, partition_budget - snapshot_s2_pair_incumbent);
                if (guaranteed_k1_min <= guaranteed_k1_max) {
                    return commitSplitDecision(true, guaranteed_k1_min);
                }
            }
            if (snapshot_s1_empty && snapshot_s2_empty &&
                snapshot_s1_pair_incumbent < std::numeric_limits<int>::max() &&
                snapshot_s2_pair_incumbent < std::numeric_limits<int>::max()) {
                const int guaranteed_k1_min = std::max(partition_min_k1, snapshot_s1_pair_incumbent);
                const int guaranteed_k1_max = std::min(partition_max_k1,
                                                       partition_budget - snapshot_s2_pair_incumbent);
                if (guaranteed_k1_min <= guaranteed_k1_max) {
                    return commitSplitDecision(true, guaranteed_k1_min);
                }
            }
            if (side1_exact_lb >= 0) {
                partition_min_k1 = std::max(partition_min_k1, side1_exact_lb);
            }
            if (partition_min_k1 > partition_max_k1) {
                return commitSplitDecision(false);
            }

            const std::size_t min_budget_problem_size =
                g_active_problem_size > 0 ? g_active_problem_size : T_work.original_nodes.size();
            if (tdsPartitionMinBudgetSearchEnabledForSize(min_budget_problem_size)) {
                const int INF_BUDGET = std::numeric_limits<int>::max() / 4;
                const auto bounds_max_fail = [](const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                              const Key128& key) -> int {
                    auto it = bounds_map.find(key);
                    if (it == bounds_map.end()) return -1;
                    return it->second.max_fail;
                };
                const auto pair_max_fail = [](const Key128& key) -> int {
                    auto it = g_kbounds.find(key);
                    if (it == g_kbounds.end()) return -1;
                    return it->second.max_fail;
                };
                const auto exact_min_success = [](const std::unordered_map<Key128, KBounds, Key128Hash>& bounds_map,
                                                 const Key128& key) -> int {
                    auto it = bounds_map.find(key);
                    if (it == bounds_map.end()) return -1;
                    if (it->second.min_success < 0) return -1;
                    if (it->second.max_fail == it->second.min_success - 1) return it->second.min_success;
                    return -1;
                };
                const auto exact_min_success_both = [&](const Key128& pair_key,
                                                       const Key128& bounds_key) -> int {
                    const int state_exact = exact_min_success(g_tds_bounds, bounds_key);
                    const int pair_exact = exact_min_success(g_kbounds, pair_key);
                    return std::max(state_exact, pair_exact);
                };
                auto makeSideRangeCacheKey = [](const Key128& side_pair_key,
                                               const Key128& side_bounds_key,
                                               const int side_min_budget,
                                               const int side_max_budget) {
                    return FreeEdgePartitionSideBudgetRangeCacheKey{
                        side_pair_key,
                        side_bounds_key,
                        static_cast<std::uint32_t>(side_min_budget),
                        static_cast<std::uint32_t>(side_max_budget)};
                };

                const auto minimalFeasibleBudget = [&](const VectorRangeTreeMap& side_bar,
                                                      const VectorRangeTreeMap& side_end,
                                                      const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& side_pairs,
                                                      const Key128 side_bounds_key,
                                                      const Key128 side_pair_key,
                                                      const int side_conflicts,
                                                      const int side_pair_incumbent,
                                                      const bool side_empty,
                                                      int side_min_budget,
                                                      int side_max_budget,
                                                      const bool is_side1) -> int {
                    if (side_min_budget > side_max_budget) {
                        return INF_BUDGET;
                    }
                    side_min_budget = std::max(side_min_budget, side_conflicts);
                    const int side_max_fail = std::max(bounds_max_fail(g_tds_bounds, side_bounds_key),
                                                      pair_max_fail(side_pair_key));
                    if (side_max_fail >= 0) {
                        side_min_budget = std::max(side_min_budget, side_max_fail + 1);
                    }
                    const int side_exact_lb = exact_min_success_both(side_pair_key, side_bounds_key);
                    if (side_exact_lb >= 0) {
                        if (side_exact_lb > side_max_budget) {
                            return INF_BUDGET;
                        }
                        side_min_budget = std::max(side_min_budget, side_exact_lb);
                    }

                    if (side_min_budget > side_max_budget) {
                        return INF_BUDGET;
                    }
                    if (side_empty && side_pair_incumbent < std::numeric_limits<int>::max() &&
                        side_pair_incumbent <= side_max_budget) {
                        if (side_pair_incumbent <= side_min_budget) {
                            return side_min_budget;
                        }
                        side_max_budget = std::min(side_max_budget, side_pair_incumbent);
                    }

                    const auto state_cache_key = makeSideBudgetCacheKey(side_pair_key, side_bounds_key);
                    const auto state_cache_it = g_free_edge_partition_side_budget_cache.find(state_cache_key);
                    if (state_cache_it != g_free_edge_partition_side_budget_cache.end() &&
                        state_cache_it->second.computed) {
                        auto& state_cache_entry = state_cache_it->second;
                        const int state_not_feasible = state_cache_entry.known_not_feasible_below;
                        if (state_not_feasible >= side_max_budget) {
                            if (g_profile.enabled) {
                                g_profile.tds_partition_side_budget_cache_hits++;
                            }
                            return INF_BUDGET;
                        }
                        if (state_not_feasible + 1 > side_min_budget) {
                            side_min_budget = state_not_feasible + 1;
                        }
                        const int known_feasible = state_cache_entry.known_min_feasible_budget;
                        if (known_feasible < std::numeric_limits<int>::max()) {
                            if (known_feasible <= side_min_budget) {
                                if (g_profile.enabled) {
                                    g_profile.tds_partition_side_budget_cache_hits++;
                                }
                                return side_min_budget;
                            }
                            if (known_feasible <= side_max_budget) {
                                if (state_not_feasible >= known_feasible - 1) {
                                    if (g_profile.enabled) {
                                        g_profile.tds_partition_side_budget_cache_hits++;
                                    }
                                    return known_feasible;
                                }
                                side_max_budget = std::min(side_max_budget, known_feasible);
                            }
                        }
                        if (side_min_budget > side_max_budget) {
                            return INF_BUDGET;
                        }
                    }

                    FreeEdgePartitionSideBudgetRangeCacheKey cache_key =
                        makeSideRangeCacheKey(side_pair_key, side_bounds_key, side_min_budget, side_max_budget);
                    const auto cached_side_range_it = g_free_edge_partition_side_budget_range_cache.find(cache_key);
                    if (cached_side_range_it != g_free_edge_partition_side_budget_range_cache.end() &&
                        cached_side_range_it->second.computed) {
                        if (g_profile.enabled) {
                            g_profile.tds_partition_side_budget_cache_hits++;
                        }
                        auto& state_cache_entry = g_free_edge_partition_side_budget_cache
                                                     .try_emplace(state_cache_key,
                                                                  FreeEdgePartitionSideBudgetCacheEntry{})
                                                     .first->second;
                        state_cache_entry.computed = true;
                        const int cached_min_feasible = cached_side_range_it->second.min_feasible_budget;
                        if (cached_min_feasible == INF_BUDGET) {
                            state_cache_entry.known_not_feasible_below = std::max(state_cache_entry.known_not_feasible_below,
                                                                                  side_max_budget);
                        } else {
                            state_cache_entry.known_min_feasible_budget =
                                std::min(state_cache_entry.known_min_feasible_budget, cached_min_feasible);
                            state_cache_entry.known_not_feasible_below = std::max(state_cache_entry.known_not_feasible_below,
                                                                                  cached_min_feasible - 1);
                        }
                        return cached_side_range_it->second.min_feasible_budget;
                    }
                    if (g_profile.enabled) {
                        g_profile.tds_partition_side_budget_cache_misses++;
                    }

                    auto isFeasible = [&](const int budget) {
                        return checkSidePartitionFeasibility(side_bar,
                                                            side_end,
                                                            side_pairs,
                                                            side_bounds_key,
                                                            side_pair_key,
                                                            side_pair_incumbent,
                                                            side_empty,
                                                            budget,
                                                            is_side1);
                    };

                    int lo = side_min_budget;
                    int hi = side_max_budget;
                    if (lo > hi) {
                        return INF_BUDGET;
                    }
                    int best = INF_BUDGET;
                    int last_infeasible = side_min_budget - 1;
                    if (tdsPartitionLowerFirstEnabled()) {
                        const bool feasible = isFeasible(lo);
                        if (g_profile.enabled) {
                            g_profile.partition_budget_iterations++;
                        }
                        if (feasible) {
                            best = lo;
                            hi = lo - 1;
                        } else {
                            last_infeasible = lo;
                            ++lo;
                        }
                    }
                    while (lo <= hi) {
                        const int mid = lo + ((hi - lo) >> 1);
                        const bool feasible = isFeasible(mid);
                        if (g_profile.enabled) {
                            g_profile.partition_budget_iterations++;
                        }
                        if (feasible) {
                            best = std::min(best, mid);
                            hi = mid - 1;
                        } else {
                            if (mid > last_infeasible) {
                                last_infeasible = mid;
                            }
                            lo = mid + 1;
                        }
                    }

                    auto& cache_entry = g_free_edge_partition_side_budget_range_cache
                                            .try_emplace(cache_key, FreeEdgePartitionSideBudgetRangeCacheEntry{})
                                            .first->second;
                    cache_entry.computed = true;
                    cache_entry.min_feasible_budget = best;

                    auto& state_cache_entry = g_free_edge_partition_side_budget_cache
                                                 .try_emplace(state_cache_key,
                                                              FreeEdgePartitionSideBudgetCacheEntry{})
                                                 .first->second;
                    state_cache_entry.computed = true;
                    if (best == INF_BUDGET) {
                        state_cache_entry.known_not_feasible_below = std::max(state_cache_entry.known_not_feasible_below,
                                                                              last_infeasible);
                    } else {
                        state_cache_entry.known_min_feasible_budget = std::min(state_cache_entry.known_min_feasible_budget,
                                                                              best);
                        state_cache_entry.known_not_feasible_below = std::max(state_cache_entry.known_not_feasible_below,
                                                                              best - 1);
                    }
                    return best;
                };

                ScopedTimer _tloop(g_profile.enabled ? &g_profile.time_partition_budget_loop_ms : nullptr);
                const int min_side_mode = tdsPartitionMinSideMode();
                const bool search_side2_first =
                    (min_side_mode == 2) ||
                    (min_side_mode == 3 && side2_lb > side1_lb) ||
                    (min_side_mode == 4 && snapshot_conf2 > snapshot_conf1) ||
                    (min_side_mode == 5 && snapshot_n2 > snapshot_n1) ||
                    (min_side_mode == 6 && snapshot_S2.size() > snapshot_S1.size());
                if (search_side2_first) {
                    const int min_side2_budget = minimalFeasibleBudget(snapshot_bar2,
                                                                       snapshot_end2,
                                                                       snapshot_S2,
                                                                       snapshot_s2_bounds_key,
                                                                       snapshot_s2_pair_key,
                                                                       snapshot_conf2,
                                                                       snapshot_s2_pair_incumbent,
                                                                       snapshot_s2_empty,
                                                                       partition_budget - partition_max_k1,
                                                                       partition_budget - partition_min_k1,
                                                                       false);
                    if (min_side2_budget == INF_BUDGET) {
                        return commitSplitDecision(false);
                    }
                    const int remaining_after_side2 = partition_budget - min_side2_budget;
                    if (remaining_after_side2 < side1_lb) {
                        return commitSplitDecision(false);
                    }
                    if (side1_exact_lb >= 0) {
                        return commitSplitDecision(side1_exact_lb <= remaining_after_side2,
                                                   remaining_after_side2);
                    }
                    if (snapshot_s1_empty &&
                        snapshot_s1_pair_incumbent < std::numeric_limits<int>::max() &&
                        snapshot_s1_pair_incumbent <= remaining_after_side2) {
                        return commitSplitDecision(true, remaining_after_side2);
                    }

                    const bool side1_ok = checkSidePartitionFeasibility(snapshot_bar1,
                                                                        snapshot_end1,
                                                                        snapshot_S1,
                                                                        snapshot_s1_bounds_key,
                                                                        snapshot_s1_pair_key,
                                                                        snapshot_s1_pair_incumbent,
                                                                        snapshot_s1_empty,
                                                                        remaining_after_side2,
                                                                        true);
                    return commitSplitDecision(side1_ok, side1_ok ? remaining_after_side2 : -1);
                }

                const int min_side1_budget = minimalFeasibleBudget(snapshot_bar1,
                                                                   snapshot_end1,
                                                                   snapshot_S1,
                                                                   snapshot_s1_bounds_key,
                                                                   snapshot_s1_pair_key,
                                                                   snapshot_conf1,
                                                                   snapshot_s1_pair_incumbent,
                                                                   snapshot_s1_empty,
                                                                   partition_min_k1,
                                                                   partition_max_k1,
                                                                   true);
                if (min_side1_budget == INF_BUDGET) {
                    return commitSplitDecision(false);
                }
                const int remaining_after_side1 = partition_budget - min_side1_budget;
                if (remaining_after_side1 < side2_lb) {
                    return commitSplitDecision(false);
                }
                if (side2_exact_lb >= 0) {
                    return commitSplitDecision(side2_exact_lb <= remaining_after_side1, min_side1_budget);
                }
                if (snapshot_s2_empty &&
                    snapshot_s2_pair_incumbent < std::numeric_limits<int>::max() &&
                    snapshot_s2_pair_incumbent <= remaining_after_side1) {
                    return commitSplitDecision(true, min_side1_budget);
                }

                const bool side2_ok = checkSidePartitionFeasibility(snapshot_bar2,
                                                                    snapshot_end2,
                                                                    snapshot_S2,
                                                                    snapshot_s2_bounds_key,
                                                                    snapshot_s2_pair_key,
                                                                    snapshot_s2_pair_incumbent,
                                                                    snapshot_s2_empty,
                                                                    remaining_after_side1,
                                                                    false);
                return commitSplitDecision(side2_ok, side2_ok ? min_side1_budget : -1);
            }

            int dynamic_side1_fail = sideKnownFailMax(snapshot_s1_pair_key, snapshot_s1_bounds_key);
            int dynamic_side2_fail = sideKnownFailMax(snapshot_s2_pair_key, snapshot_s2_bounds_key);
            int dynamic_side1_exact_lb = exactMinSuccessOrZero(snapshot_s1_pair_key, snapshot_s1_bounds_key);
            int dynamic_side2_exact_lb = exactMinSuccessOrZero(snapshot_s2_pair_key, snapshot_s2_bounds_key);

            const auto recomputePartitionWindow = [&]() -> bool {
                partition_min_k1 = std::max(0, side1_lb);
                if (dynamic_side1_fail >= 0) {
                    partition_min_k1 = std::max(partition_min_k1, dynamic_side1_fail + 1);
                }
                if (dynamic_side1_exact_lb >= 0) {
                    partition_min_k1 = std::max(partition_min_k1, dynamic_side1_exact_lb);
                }
                partition_max_k1 = std::min(partition_budget, partition_budget - side2_lb);
                if (dynamic_side2_fail >= 0) {
                    partition_max_k1 = std::min(partition_max_k1, partition_budget - (dynamic_side2_fail + 1));
                }
                if (dynamic_side2_exact_lb >= 0) {
                    partition_max_k1 = std::min(partition_max_k1, partition_budget - dynamic_side2_exact_lb);
                }
                return partition_min_k1 <= partition_max_k1;
            };

            const auto advanceK1AfterFail = [](const int current_k1, const int next_min_k1,
                                              const int fail_limit, const int exact_lb) {
                int next_k1 = current_k1 + 1;
                if (next_k1 < next_min_k1) {
                    next_k1 = next_min_k1;
                }
                if (fail_limit >= 0 && fail_limit + 1 > next_k1) {
                    next_k1 = fail_limit + 1;
                }
                if (exact_lb >= 0 && exact_lb > next_k1) {
                    next_k1 = exact_lb;
                }
                return next_k1;
            };

            if (!recomputePartitionWindow()) {
                return commitSplitDecision(false);
            }
            if (side1_exact_lb >= 0 && side2_exact_lb >= 0) {
                const int guaranteed_k1_min = std::max(partition_min_k1, side1_exact_lb);
                const int guaranteed_k1_max = std::min(partition_max_k1, partition_budget - side2_exact_lb);
                if (guaranteed_k1_min <= guaranteed_k1_max) {
                    return commitSplitDecision(true, guaranteed_k1_min);
                }
            }

            {
                ScopedTimer _tloop(g_profile.enabled ? &g_profile.time_partition_budget_loop_ms : nullptr);
                for (int k1 = partition_min_k1; k1 <= partition_max_k1;) {
                    if (g_profile.enabled) g_profile.partition_budget_iterations++;
                    if (!recomputePartitionWindow()) {
                        return commitSplitDecision(false);
                    }
                    if (k1 < partition_min_k1) {
                        k1 = partition_min_k1;
                    }
                    const int k2 = partition_budget - k1;
                    if (k2 < 0) {
                        return commitSplitDecision(false);
                    }
                    if (dynamic_side2_fail >= 0 && k2 <= dynamic_side2_fail) {
                        return commitSplitDecision(false);
                    }
                    if (dynamic_side1_fail >= 0 && k1 <= dynamic_side1_fail) {
                        k1 = std::max(partition_min_k1, dynamic_side1_fail + 1);
                        continue;
                    }
                    if (dynamic_side1_exact_lb >= 0 && k1 < dynamic_side1_exact_lb) {
                        k1 = std::max(partition_min_k1, dynamic_side1_exact_lb);
                        continue;
                    }
                    if (dynamic_side2_exact_lb >= 0 && k2 < dynamic_side2_exact_lb) {
                        return commitSplitDecision(false);
                    }

                    const bool side1_guaranteed = (dynamic_side1_exact_lb >= 0 && k1 >= dynamic_side1_exact_lb);
                    const bool side2_guaranteed = (dynamic_side2_exact_lb >= 0 && k2 >= dynamic_side2_exact_lb);
                    if (side1_guaranteed && side2_guaranteed) {
                        return commitSplitDecision(true, k1);
                    }

                    if (side1_guaranteed) {
                        const bool s2_ok = checkSidePartitionFeasibility(snapshot_bar2,
                                                                        snapshot_end2,
                                                                        snapshot_S2,
                                                                        snapshot_s2_bounds_key,
                                                                        snapshot_s2_pair_key,
                                                                        snapshot_s2_pair_incumbent,
                                                                        snapshot_s2_empty,
                                                                        k2,
                                                                        false);
                        if (!s2_ok) {
                            const int s2_fail = sideKnownFailMax(snapshot_s2_pair_key, snapshot_s2_bounds_key);
                            dynamic_side2_fail = std::max(dynamic_side2_fail, s2_fail);
                            dynamic_side2_exact_lb = std::max(dynamic_side2_exact_lb,
                                                              exactMinSuccessOrZero(snapshot_s2_pair_key,
                                                                               snapshot_s2_bounds_key));
                            if (dynamic_side2_fail >= k2 || (dynamic_side2_exact_lb >= 0 && k2 < dynamic_side2_exact_lb)) {
                                return commitSplitDecision(false);
                            }
                            if (!recomputePartitionWindow()) {
                                return commitSplitDecision(false);
                            }
                            k1 = advanceK1AfterFail(k1,
                                                    partition_min_k1,
                                                    dynamic_side2_fail,
                                                    dynamic_side2_exact_lb);
                            continue;
                        }
                        dynamic_side2_exact_lb = std::max(dynamic_side2_exact_lb,
                                                          exactMinSuccessOrZero(snapshot_s2_pair_key,
                                                                               snapshot_s2_bounds_key));
                        return commitSplitDecision(true, k1);
                    }

                    if (side2_guaranteed) {
                        const bool s1_ok = checkSidePartitionFeasibility(snapshot_bar1,
                                                                        snapshot_end1,
                                                                        snapshot_S1,
                                                                        snapshot_s1_bounds_key,
                                                                        snapshot_s1_pair_key,
                                                                        snapshot_s1_pair_incumbent,
                                                                        snapshot_s1_empty,
                                                                        k1,
                                                                        true);
                        if (!s1_ok) {
                            const int s1_fail = sideKnownFailMax(snapshot_s1_pair_key, snapshot_s1_bounds_key);
                            dynamic_side1_fail = std::max(dynamic_side1_fail, s1_fail);
                            dynamic_side1_exact_lb = std::max(dynamic_side1_exact_lb,
                                                              exactMinSuccessOrZero(snapshot_s1_pair_key,
                                                                                   snapshot_s1_bounds_key));
                            if (!recomputePartitionWindow()) {
                                return commitSplitDecision(false);
                            }
                            k1 = advanceK1AfterFail(k1,
                                                    partition_min_k1,
                                                    dynamic_side1_fail,
                                                    dynamic_side1_exact_lb);
                            continue;
                        }
                        dynamic_side1_exact_lb = std::max(dynamic_side1_exact_lb,
                                                          exactMinSuccessOrZero(snapshot_s1_pair_key,
                                                                               snapshot_s1_bounds_key));
                        return commitSplitDecision(true, k1);
                    }

                    const int side1_min_for_split = std::max({side1_lb, dynamic_side1_fail + 1, dynamic_side1_exact_lb});
                    const int side2_min_for_split = std::max({side2_lb, dynamic_side2_fail + 1, dynamic_side2_exact_lb});
                    const bool check_side1_first = (k1 - side1_min_for_split) <= (k2 - side2_min_for_split);

                    if (check_side1_first) {
                        const bool s1_ok = checkSidePartitionFeasibility(snapshot_bar1,
                                                                        snapshot_end1,
                                                                        snapshot_S1,
                                                                        snapshot_s1_bounds_key,
                                                                        snapshot_s1_pair_key,
                                                                        snapshot_s1_pair_incumbent,
                                                                        snapshot_s1_empty,
                                                                        k1,
                                                                        true);
                        if (!s1_ok) {
                            const int s1_fail = sideKnownFailMax(snapshot_s1_pair_key, snapshot_s1_bounds_key);
                            dynamic_side1_fail = std::max(dynamic_side1_fail, s1_fail);
                            dynamic_side1_exact_lb = std::max(dynamic_side1_exact_lb,
                                                              exactMinSuccessOrZero(snapshot_s1_pair_key,
                                                                               snapshot_s1_bounds_key));
                            if (!recomputePartitionWindow()) {
                                return commitSplitDecision(false);
                            }
                            k1 = advanceK1AfterFail(k1,
                                                    partition_min_k1,
                                                    dynamic_side1_fail,
                                                    dynamic_side1_exact_lb);
                            continue;
                        }
                        dynamic_side1_exact_lb = std::max(dynamic_side1_exact_lb,
                                                          exactMinSuccessOrZero(snapshot_s1_pair_key,
                                                                               snapshot_s1_bounds_key));

                        const bool s2_ok = checkSidePartitionFeasibility(snapshot_bar2,
                                                                        snapshot_end2,
                                                                        snapshot_S2,
                                                                        snapshot_s2_bounds_key,
                                                                        snapshot_s2_pair_key,
                                                                        snapshot_s2_pair_incumbent,
                                                                        snapshot_s2_empty,
                                                                        k2,
                                                                        false);
                        if (!s2_ok) {
                            const int s2_fail = sideKnownFailMax(snapshot_s2_pair_key, snapshot_s2_bounds_key);
                            dynamic_side2_fail = std::max(dynamic_side2_fail, s2_fail);
                            dynamic_side2_exact_lb = std::max(dynamic_side2_exact_lb,
                                                              exactMinSuccessOrZero(snapshot_s2_pair_key,
                                                                               snapshot_s2_bounds_key));
                            if (dynamic_side2_fail >= k2 || (dynamic_side2_exact_lb >= 0 && k2 < dynamic_side2_exact_lb)) {
                                return commitSplitDecision(false);
                            }
                            if (!recomputePartitionWindow()) {
                                return commitSplitDecision(false);
                            }
                            k1 = advanceK1AfterFail(k1,
                                                    partition_min_k1,
                                                    dynamic_side2_fail,
                                                    dynamic_side2_exact_lb);
                            continue;
                        }
                        dynamic_side2_exact_lb = std::max(dynamic_side2_exact_lb,
                                                          exactMinSuccessOrZero(snapshot_s2_pair_key,
                                                                               snapshot_s2_bounds_key));
                        return commitSplitDecision(true, k1);
                    }

                    const bool s2_ok = checkSidePartitionFeasibility(snapshot_bar2,
                                                                    snapshot_end2,
                                                                    snapshot_S2,
                                                                    snapshot_s2_bounds_key,
                                                                    snapshot_s2_pair_key,
                                                                    snapshot_s2_pair_incumbent,
                                                                    snapshot_s2_empty,
                                                                    k2,
                                                                    false);
                    if (!s2_ok) {
                        const int s2_fail = sideKnownFailMax(snapshot_s2_pair_key, snapshot_s2_bounds_key);
                        dynamic_side2_fail = std::max(dynamic_side2_fail, s2_fail);
                        dynamic_side2_exact_lb = std::max(dynamic_side2_exact_lb,
                                                          exactMinSuccessOrZero(snapshot_s2_pair_key,
                                                                               snapshot_s2_bounds_key));
                        if (dynamic_side2_fail >= k2) {
                            return commitSplitDecision(false);
                        }
                        if (!recomputePartitionWindow()) {
                            return commitSplitDecision(false);
                        }
                        k1 = advanceK1AfterFail(k1,
                                                partition_min_k1,
                                                dynamic_side2_fail,
                                                dynamic_side2_exact_lb);
                        continue;
                    }

                    const bool s1_ok = checkSidePartitionFeasibility(snapshot_bar1,
                                                                    snapshot_end1,
                                                                    snapshot_S1,
                                                                    snapshot_s1_bounds_key,
                                                                    snapshot_s1_pair_key,
                                                                    snapshot_s1_pair_incumbent,
                                                                    snapshot_s1_empty,
                                                                    k1,
                                                                    true);
                if (!s1_ok) {
                    const int s1_fail = sideKnownFailMax(snapshot_s1_pair_key, snapshot_s1_bounds_key);
                    dynamic_side1_fail = std::max(dynamic_side1_fail, s1_fail);
                    dynamic_side1_exact_lb = std::max(dynamic_side1_exact_lb,
                                                          exactMinSuccessOrZero(snapshot_s1_pair_key,
                                                                               snapshot_s1_bounds_key));
                    if (!recomputePartitionWindow()) {
                        return commitSplitDecision(false);
                    }
                    k1 = advanceK1AfterFail(k1,
                                            partition_min_k1,
                                            dynamic_side1_fail,
                                            dynamic_side1_exact_lb);
                    continue;
                }
                    dynamic_side1_exact_lb = std::max(dynamic_side1_exact_lb,
                                                      exactMinSuccessOrZero(snapshot_s1_pair_key,
                                                                           snapshot_s1_bounds_key));
                    return commitSplitDecision(true, k1);
                }
            }
            return commitSplitDecision(false);
        } catch (const SearchAbortException&) {
            throw;
        } catch (...) {
            debugPrint("TreeDistS: Exception during free edge handling");
            return commitResult(false);
        }
    }

    //                     For each nonempty independent subset I ⊆ ⋃ S (no two edges in I share a node):"
    debugPrint("TreeDistS: No free edge, implementing S branching (step 2)");

    if (k_work <= 0) {
        debugPrint("TreeDistS: No budget left");
        return commitResult(false);
    }

    if (S_work.empty()) {
        // No constraints from S, try any rotation
        if (g_profile.enabled) {
            g_profile.s_empty_branch_calls++;
        }
        debugPrint("TreeDistS: S is empty, trying any rotation");
        if (emptySCoreReduceEnabled()) {
            EmptySCoreSummary core_summary = analyzeEmptySCore(T_work, T_end);
            if (core_summary.fully_exact) {
                return commitResult(core_summary.exact_total <= k_work);
            }
            if (core_summary.lower_bound > k_work) {
                return commitResult(false);
            }
        }
        if (emptySMotifReduceEnabled()) {
            bool motif_handled = false;
            bool motif_result = false;
            if (tryExactEmptySMotifReduction(T_work, T_end, k_work, motif_handled, motif_result) &&
                motif_handled) {
                return commitResult(motif_result);
            }
        }
        if (k_work >= 1 && tdsCommonEdgeRotatePairEnabled()) {
            const int min_common_edges = tdsCommonEdgeMinCount();
            const int common_edge_budget_cap = std::max(4, min_common_edges * 5);
            if (k_work <= common_edge_budget_cap) {
                const bool has_common_edge_budget = pairHasCommonEdgesCached(T_work, T_end, min_common_edges);
                if (has_common_edge_budget) {
                    bool handled_common = false;
                    if (tryCommonEdgeDecomposeCached(T_work, T_end, k_work, handled_common)) {
                        return commitResult(true);
                    }
                }
                if (has_common_edge_budget && k_work >= 2 &&
                    tryCommonEdgeDecomposeAfterOneRotationEachCached(T_work, T_end, k_work)) {
                    return commitResult(true);
                }
            }
        }
        const int current_conflicts = S_work.empty()
                                          ? cachedConflictCountWithKey(T_work, T_end, work_pair_key)
                                          : cachedConflictCount(T_work, T_end);
        const bool must_drop_conflict_now = (k_work == current_conflicts);
        const int child_k = k_work - 1;
        const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> empty_S;
        struct RotationCandidate {
            std::pair<int,int> edge;
            int next_conflicts = INT_MAX;
            int next_conflicts_bucket = INT_MAX;
            int free_hint = -1;
            int child_lb = INT_MAX;
            Key128 child_exact_key{};
            Key128 child_bounds_key{};
            Key128 child_pair_key{};
            std::uint64_t partition_signature = 0;
        };
        const bool use_signature_dedupe = emptySDedupEnabled();
        const int order_mode = emptySOrderMode();
        const bool use_free_hint = (g_empty_s_disable_free_hint_depth == 0);
        const bool use_incumbent_prune = incumbentPruneEnabled();
        const auto computeFreeHint = [&](RotationCandidate& candidate, const VectorRangeTreeMap& probe_state) {
            auto [probe_has_free, probe_free_edge] =
                cachedFindFreeEdge(probe_state, T_end, candidate.child_pair_key);
            candidate.free_hint = probe_has_free ? 1 : 0;
            if (candidate.free_hint && probe_free_edge.first >= 0 && probe_free_edge.second >= 0) {
                int a = std::min(probe_free_edge.first, probe_free_edge.second);
                int b = std::max(probe_free_edge.first, probe_free_edge.second);
                candidate.partition_signature =
                    (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
                    static_cast<std::uint64_t>(static_cast<std::uint32_t>(b));
            } else {
                candidate.partition_signature = 0;
            }
        };
        const auto candidateDominates = [&](RotationCandidate& a, const RotationCandidate& b) {
            if (a.child_lb != b.child_lb) {
                return a.child_lb < b.child_lb;
            }
            if (a.next_conflicts != b.next_conflicts) {
                return a.next_conflicts < b.next_conflicts;
            }
            if (a.free_hint != b.free_hint) {
                return a.free_hint > b.free_hint;
            }
            return false;
        };
        const auto candidateStateDominates = [&](RotationCandidate& a, const RotationCandidate& b) {
            if (a.child_lb != b.child_lb) {
                return a.child_lb < b.child_lb;
            }
            if (a.next_conflicts != b.next_conflicts) {
                return a.next_conflicts < b.next_conflicts;
            }
            if (a.free_hint != b.free_hint) {
                return a.free_hint > b.free_hint;
            }
            return false;
        };
        struct EmptySChildStateSignature {
            Key128 child_pair_key{};
            Key128 child_bounds_key{};
            int next_conflicts_bucket = INT_MAX;
            bool operator==(const EmptySChildStateSignature& o) const {
                return child_pair_key == o.child_pair_key &&
                       child_bounds_key == o.child_bounds_key &&
                       next_conflicts_bucket == o.next_conflicts_bucket;
            }
        };
        struct EmptySChildStateSignatureHash {
            std::size_t operator()(const EmptySChildStateSignature& s) const {
                std::size_t h = Key128Hash{}(s.child_pair_key);
                h ^= Key128Hash{}(s.child_bounds_key) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                h ^= std::hash<int>{}(s.next_conflicts_bucket) + 0x517cc1b727220a95ULL + (h << 5) + (h >> 3);
                return h;
            }
        };
        const auto populateCandidate = [&](RotationCandidate& cand,
                                           const VectorRangeTreeMap& candidate_tree,
                                           const int next_conflicts_hint) -> bool {
            cand.next_conflicts = (next_conflicts_hint >= 0)
                                      ? next_conflicts_hint
                                      : cachedConflictCount(candidate_tree, T_end);
            if (cand.next_conflicts < 0) {
                return false;
            }
            cand.child_lb = std::max(requiredBudgetFromBounds(g_tds_bounds, cand.child_bounds_key),
                                   requiredBudgetFromBounds(g_kbounds, cand.child_pair_key));
            cand.child_lb = std::max(cand.next_conflicts, cand.child_lb);
            return true;
        };
        auto edges = getInternalEdges(T_work);
        const auto& s_empty_target_edges = cachedTargetEdgeSet(T_end);
        std::vector<RotationCandidate> candidates;
        candidates.reserve(edges.size());
        std::vector<RotationCandidate> plateau_candidates;
        plateau_candidates.reserve(edges.size());
        auto byPromise = [order_mode](const RotationCandidate& a, const RotationCandidate& b) {
            if (order_mode == 1 && a.free_hint != b.free_hint) {
                return a.free_hint > b.free_hint;
            }
            if (order_mode == 2 && a.next_conflicts != b.next_conflicts) {
                return a.next_conflicts < b.next_conflicts;
            }
            if (a.child_lb != b.child_lb) {
                return a.child_lb < b.child_lb;
            }
            if (a.next_conflicts != b.next_conflicts) {
                return a.next_conflicts < b.next_conflicts;
            }
            if (a.free_hint != b.free_hint) {
                return a.free_hint > b.free_hint;
            }
            if (a.edge.first != b.edge.first) {
                return a.edge.first < b.edge.first;
            }
            return a.edge.second < b.edge.second;
        };
        auto rotateEdgeState = [&](int parent, int child, int& rotate_dir) -> bool {
            if (T_work.getLeftChild(parent) == child) {
                T_work.rotateRight(parent);
                rotate_dir = 1;
                return true;
            }
            if (T_work.getRightChild(parent) == child) {
                T_work.rotateLeft(parent);
                rotate_dir = -1;
                return true;
            }
            rotate_dir = 0;
            return false;
        };
        auto undoRotateState = [&](int child, int rotate_dir) {
            if (rotate_dir == 1) {
                T_work.rotateLeft(child);
            } else if (rotate_dir == -1) {
                T_work.rotateRight(child);
            }
        };
        const auto refreshDeferredFreeHint = [&](RotationCandidate& cand,
                                                EmptySChildCacheItem* cache_item) -> bool {
            if (cand.free_hint >= 0) {
                return true;
            }
            if (!use_free_hint) {
                cand.free_hint = 0;
                cand.partition_signature = 0;
                if (cache_item) {
                    cache_item->free_hint = cand.free_hint;
                    cache_item->partition_signature = cand.partition_signature;
                }
                return true;
            }
            int rotate_dir = 0;
            bool rotated = false;
            const int parent = cand.edge.first;
            const int child = cand.edge.second;
            auto rollbackRotation = [&]() {
                if (rotated) {
                    undoRotateState(child, rotate_dir);
                    rotated = false;
                }
            };
            try {
                rotated = rotateEdgeState(parent, child, rotate_dir);
                if (!rotated) {
                    rollbackRotation();
                    return false;
                }
                computeFreeHint(cand, T_work);
                rollbackRotation();
                if (cache_item) {
                    cache_item->free_hint = cand.free_hint;
                    cache_item->partition_signature = cand.partition_signature;
                }
                return true;
            } catch (const SearchAbortException&) {
                throw;
            } catch (...) {
                rollbackRotation();
                return false;
            }
        };
        if (use_signature_dedupe) {
            const EmptySChildCacheKey cache_key{work_pair_key};
            std::vector<RotationCandidate> cached_candidates;
            bool early_commit = false;
            auto buildCacheAndCollect = [&](bool build_new_cache) -> std::vector<RotationCandidate> {
                std::vector<std::pair<EmptySChildStateSignature, RotationCandidate>> best_signature_states;
                best_signature_states.reserve(edges.size());
                std::vector<Key128> seen_child_pair_keys;
                seen_child_pair_keys.reserve(edges.size());
	                for (const auto& edge : edges) {
	                    int rotate_dir = 0;
	                    bool rotated = false;
	                    int parent = edge.first;
	                    int child = edge.second;
	                    std::array<RP, 4> old_rotation_edges{};
	                    int old_rotation_edge_count = 0;
	                    auto rollbackRotation = [&]() {
	                        if (rotated) {
	                            undoRotateState(child, rotate_dir);
                            rotated = false;
                        }
	                    };
	                    try {
	                        old_rotation_edge_count =
	                            collectOutgoingRangeEdgesForNodes(T_work, parent, child, old_rotation_edges);
	                        rotated = rotateEdgeState(parent, child, rotate_dir);
	                        if (!rotated) {
	                            rollbackRotation();
                            continue;
                        }

	                        const int next_conflicts_hint =
	                            rotatedConflictCountFromLocalDelta(current_conflicts,
	                                                               old_rotation_edges,
	                                                               old_rotation_edge_count,
	                                                               T_work,
	                                                               parent,
	                                                               child,
	                                                               s_empty_target_edges);
	                        if (next_conflicts_hint > child_k) {
	                            rollbackRotation();
	                            continue;
	                        }
	                        if (!build_new_cache && must_drop_conflict_now &&
	                            next_conflicts_hint >= current_conflicts) {
	                            rollbackRotation();
	                            if (g_profile.enabled) {
	                                g_profile.s_empty_must_drop_rejects++;
	                            }
	                            continue;
	                        }

	                        EmptySKeys child_keys = makeEmptySKeys(T_work, T_end, child_k);
                        Key128 child_pair_key = child_keys.pair_key;
                        const bool repeated_pair_state =
                            std::find(seen_child_pair_keys.begin(),
                                      seen_child_pair_keys.end(),
                                      child_pair_key) != seen_child_pair_keys.end();
                        if (repeated_pair_state) {
                            rollbackRotation();
                            if (g_profile.enabled) {
                                g_profile.s_empty_duplicate_child_states++;
                            }
                            continue;
                        }
                        seen_child_pair_keys.push_back(child_pair_key);

                        const int child_pair_incumbent = use_incumbent_prune
                                                            ? pairIncumbentUpper(child_pair_key)
                                                            : std::numeric_limits<int>::max();
                        if (child_pair_incumbent <= child_k) {
                            rollbackRotation();
                            if (g_profile.enabled) {
                                g_profile.s_empty_incumbent_rejects++;
                            }
                            early_commit = true;
                            break;
                        }

                        Key128 child_bounds_key = child_keys.bounds_key;
                        Key128 child_exact_key = child_keys.exact_key;

                        auto child_memo_it = g_memo_tds.find(child_exact_key);
                        if (child_memo_it != g_memo_tds.end()) {
                            if (child_memo_it->second) {
                                rollbackRotation();
                                early_commit = true;
                                break;
                            }
                        }
                        bool child_known_tds = false;
                        if (tryBoundsPrune(g_tds_bounds, child_bounds_key, child_k, child_known_tds)) {
                            if (g_profile.enabled) {
                                g_profile.bounds_hits_tds++;
                            }
                            if (child_known_tds) {
                                rollbackRotation();
                                early_commit = true;
                                break;
                            }
                        }

                        RotationCandidate cand{edge,
                                              INT_MAX,
                                              INT_MAX,
                                              -1,
                                              INT_MAX,
                                              child_exact_key,
	                                              child_bounds_key,
	                                              child_pair_key,
	                                              0};
		                        if (!populateCandidate(cand, T_work, next_conflicts_hint)) {
		                            rollbackRotation();
		                            continue;
		                        }
                        if (!build_new_cache && cand.child_lb > child_k) {
                            rollbackRotation();
                            continue;
                        }
                        const bool child_within_budget = cand.child_lb <= child_k;
                        if (cand.next_conflicts < current_conflicts) {
                            cand.free_hint = 0;
                            cand.partition_signature = 0;
                        } else if (must_drop_conflict_now) {
                            cand.free_hint = -1;
                            cand.partition_signature = 0;
                        } else if (!child_within_budget) {
                            cand.free_hint = -1;
                            cand.partition_signature = 0;
                        } else if (use_free_hint) {
                            computeFreeHint(cand, T_work);
                        } else {
                            cand.free_hint = 0;
                            cand.partition_signature = 0;
                        }
                        cand.next_conflicts_bucket = cand.next_conflicts;
                        EmptySChildStateSignature sig{
                            child_pair_key, child_bounds_key, cand.next_conflicts};
                        auto it = std::find_if(best_signature_states.begin(),
                                               best_signature_states.end(),
                                               [&](const auto& entry) {
                                                   return entry.first == sig;
                                               });
                        if (it == best_signature_states.end()) {
                            best_signature_states.emplace_back(sig, cand);
                        } else if (candidateStateDominates(cand, it->second)) {
                            it->second = cand;
                        }
                        rollbackRotation();
                    } catch (const SearchAbortException&) {
                        throw;
                    } catch (...) {
                        rollbackRotation();
                        continue;
                    }
                }
                std::vector<RotationCandidate> built_candidates;
                built_candidates.reserve(best_signature_states.size());
                if (build_new_cache) {
                    EmptySChildCacheEntry cache_entry;
                    cache_entry.items.reserve(best_signature_states.size());
                    for (auto&& entry : best_signature_states) {
                        auto& cand = entry.second;
                        cache_entry.items.push_back(EmptySChildCacheItem{
                            cand.edge,
                            cand.next_conflicts,
                            cand.next_conflicts_bucket,
                            cand.free_hint,
                            cand.child_lb,
                            child_k,
                            use_incumbent_prune
                                ? pairIncumbentUpper(cand.child_pair_key)
                                : std::numeric_limits<int>::max(),
                            cand.child_pair_key,
                            cand.child_bounds_key,
                            cand.partition_signature});
                        if (cand.child_lb <= child_k) {
                            built_candidates.push_back(cand);
                        }
                    }
                    g_empty_s_child_cache.emplace(cache_key, std::move(cache_entry));
                } else {
                    for (auto&& entry : best_signature_states) {
                        if (entry.second.child_lb <= child_k) {
                            built_candidates.push_back(std::move(entry.second));
                        }
                    }
                }
                return built_candidates;
            };
	            const bool tight_no_cache =
	                must_drop_conflict_now && emptySTightNoCacheEnabled();
	            auto cache_it = tight_no_cache ? g_empty_s_child_cache.end()
	                                           : g_empty_s_child_cache.find(cache_key);
	            if (cache_it == g_empty_s_child_cache.end()) {
	                cached_candidates = buildCacheAndCollect(!tight_no_cache);
	                if (early_commit) {
	                    return commitResult(true);
	                }
            } else {
                for (auto& item : cache_it->second.items) {
                    RotationCandidate cand;
                    cand.edge = item.edge;
                    cand.next_conflicts = item.next_conflicts;
                    cand.next_conflicts_bucket = item.next_conflicts_bucket;
                    cand.free_hint = item.free_hint;
                    cand.child_lb = item.child_lb;
                    cand.child_pair_key = item.child_pair_key;
                    cand.child_bounds_key = item.child_bounds_key;
                    cand.partition_signature = item.partition_signature;
                    if (must_drop_conflict_now && cand.next_conflicts >= current_conflicts) {
                        if (g_profile.enabled) {
                            g_profile.s_empty_must_drop_rejects++;
                        }
                        continue;
                    }
                    const int child_pair_incumbent = use_incumbent_prune
                                                        ? pairIncumbentUpper(cand.child_pair_key)
                                                        : std::numeric_limits<int>::max();
                    if (child_pair_incumbent <= child_k) {
                        return commitResult(true);
                    }
                    if (item.child_k == child_k) {
                        const int required_lb = std::max(requiredBudgetFromBounds(g_tds_bounds, cand.child_bounds_key),
                                                         requiredBudgetFromBounds(g_kbounds, cand.child_pair_key));
                        cand.child_lb = std::max({item.child_lb, required_lb, cand.next_conflicts});
                        if (cand.child_lb > child_k) {
                            continue;
                        }
                        if (!refreshDeferredFreeHint(cand, &item)) {
                            continue;
                        }
                        cached_candidates.push_back(cand);
                        continue;
                    }
                    bool child_bounds_known = false;
                    if (tryBoundsPrune(g_tds_bounds, cand.child_bounds_key, child_k, child_bounds_known)) {
                        if (g_profile.enabled) {
                            g_profile.bounds_hits_tds++;
                        }
                        if (child_bounds_known) {
                            return commitResult(true);
                        }
                        if (g_profile.enabled) {
                            g_profile.s_empty_duplicate_child_states++;
                        }
                        continue;
                    }
                    if (cand.next_conflicts < 0 || cand.next_conflicts > child_k) {
                        continue;
                    }
                    cand.child_lb = std::max(requiredBudgetFromBounds(g_tds_bounds, cand.child_bounds_key),
                                            requiredBudgetFromBounds(g_kbounds, cand.child_pair_key));
                    cand.child_lb = std::max(cand.child_lb, cand.next_conflicts);
                    if (cand.child_lb > child_k) {
                        continue;
                    }
                    if (!refreshDeferredFreeHint(cand, &item)) {
                        continue;
                    }
                    cached_candidates.push_back(cand);
                }
            }
            for (auto&& cand : cached_candidates) {
                if (must_drop_conflict_now && cand.next_conflicts >= current_conflicts) {
                    if (g_profile.enabled) {
                        g_profile.s_empty_must_drop_rejects++;
                    }
                    continue;
                }
                if (cand.next_conflicts < current_conflicts || cand.free_hint) {
                    candidates.push_back(cand);
                } else {
                    plateau_candidates.push_back(cand);
                }
            }
        } else {
            std::unordered_map<EmptySChildSignature, RotationCandidate, EmptySChildSignatureHash> best_signature_states;
            best_signature_states.reserve(edges.size() * 2 + 1);
            std::unordered_set<Key128, Key128Hash> seen_child_state_keys;
            seen_child_state_keys.reserve(edges.size());
            std::unordered_set<Key128, Key128Hash> seen_child_pair_keys;
            seen_child_pair_keys.reserve(edges.size());
	                for (const auto& edge : edges) {
	                    int rotate_dir = 0;
	                    bool rotated = false;
	                    int parent = edge.first;
	                    int child = edge.second;
	                    std::array<RP, 4> old_rotation_edges{};
	                    int old_rotation_edge_count = 0;
	                    auto rollbackRotation = [&]() {
	                        if (rotated) {
	                            undoRotateState(child, rotate_dir);
                        rotated = false;
                    }
	                };
	                try {
	                    old_rotation_edge_count =
	                        collectOutgoingRangeEdgesForNodes(T_work, parent, child, old_rotation_edges);
	                    rotated = rotateEdgeState(parent, child, rotate_dir);
	                    if (!rotated) {
	                        rollbackRotation();
                        continue;
                    }

	                    const int next_conflicts_hint =
	                        rotatedConflictCountFromLocalDelta(current_conflicts,
	                                                           old_rotation_edges,
	                                                           old_rotation_edge_count,
	                                                           T_work,
	                                                           parent,
	                                                           child,
	                                                           s_empty_target_edges);
	                    if (next_conflicts_hint > child_k) {
	                        rollbackRotation();
	                        continue;
	                    }

                    EmptySKeys child_keys = makeEmptySKeys(T_work, T_end, child_k);
                    Key128 child_pair_key = child_keys.pair_key;
                    bool unique_pair_state = seen_child_pair_keys.insert(child_pair_key).second;
                    if (!unique_pair_state) {
                        rollbackRotation();
                        if (g_profile.enabled) {
                            g_profile.s_empty_duplicate_child_states++;
                        }
                        continue;
                    }

                    const int child_pair_incumbent = use_incumbent_prune
                                                        ? pairIncumbentUpper(child_pair_key)
                                                        : std::numeric_limits<int>::max();
                    if (child_pair_incumbent <= child_k) {
                        rollbackRotation();
                        if (g_profile.enabled) {
                            g_profile.s_empty_incumbent_rejects++;
                        }
                        return commitResult(true);
                    }

                    Key128 child_bounds_key = child_keys.bounds_key;
                    Key128 child_exact_key = child_keys.exact_key;
                    bool unique_child = seen_child_state_keys.insert(child_exact_key).second;
                    if (unique_child) {
                        auto child_memo_it = g_memo_tds.find(child_exact_key);
                        if (child_memo_it != g_memo_tds.end()) {
                            if (child_memo_it->second) {
                                return commitResult(true);
                            }
                            if (g_profile.enabled) {
                                g_profile.s_empty_duplicate_child_states++;
                            }
                            rollbackRotation();
                            continue;
                        }
                        bool child_known_tds = false;
                        if (tryBoundsPrune(g_tds_bounds, child_bounds_key, child_k, child_known_tds)) {
                            if (g_profile.enabled) {
                                g_profile.bounds_hits_tds++;
                            }
                            if (child_known_tds) {
                                return commitResult(true);
                            }
                            if (g_profile.enabled) {
                                g_profile.s_empty_duplicate_child_states++;
                            }
                            rollbackRotation();
                            continue;
                        }
                        RotationCandidate cand{edge,
                                              INT_MAX,
                                              INT_MAX,
                                              -1,
                                              INT_MAX,
                                              child_exact_key,
	                                              child_bounds_key,
	                                              child_pair_key,
	                                              0};
	                        if (populateCandidate(cand, T_work, next_conflicts_hint)) {
                            if (cand.child_lb > child_k) {
                                rollbackRotation();
                                continue;
                            }
	                            if (cand.next_conflicts < current_conflicts) {
	                                cand.free_hint = 0;
	                            } else if (use_free_hint) {
	                                computeFreeHint(cand, T_work);
	                            } else {
	                                cand.free_hint = 0;
                            }
                                cand.next_conflicts_bucket = cand.next_conflicts;
                                EmptySChildSignature sig{child_pair_key,
                                                     child_bounds_key,
                                                     cand.next_conflicts_bucket};
                                auto it = best_signature_states.find(sig);
                            if (it == best_signature_states.end()) {
                                best_signature_states.emplace(sig, cand);
                            } else if (candidateDominates(cand, it->second)) {
                                it->second = cand;
                            }
                        } else {
                            unique_child = false;
                        }
                    }
                    if (!unique_child) {
                        if (g_profile.enabled) {
                            g_profile.s_empty_duplicate_child_states++;
                        }
                        rollbackRotation();
                        continue;
                    }
                    rollbackRotation();
                } catch (const SearchAbortException&) {
                    throw;
                } catch (...) {
                    rollbackRotation();
                    continue;
                }
            }

            for (auto&& it : best_signature_states) {
                auto& cand = it.second;
                if (cand.child_lb > child_k) {
                    continue;
                }
                if (cand.next_conflicts < current_conflicts || cand.free_hint) {
                    candidates.push_back(cand);
                } else {
                    plateau_candidates.push_back(cand);
                }
            }
        }

        if (g_profile.enabled) {
            g_profile.s_empty_progress_candidates += (long long)candidates.size();
            g_profile.s_empty_plateau_candidates += (long long)plateau_candidates.size();
        }
        std::sort(candidates.begin(), candidates.end(), byPromise);
        std::sort(plateau_candidates.begin(), plateau_candidates.end(), byPromise);
        candidates.insert(candidates.end(), plateau_candidates.begin(), plateau_candidates.end());

        for (const auto& cand : candidates) {
            int rotate_dir = 0;
            bool rotated = false;
            int parent = cand.edge.first;
            int child = cand.edge.second;
            auto rollbackRotation = [&]() {
                if (rotated) {
                    undoRotateState(child, rotate_dir);
                    rotated = false;
                }
            };
            try {
                rotated = rotateEdgeState(parent, child, rotate_dir);
                if (!rotated) {
                    rollbackRotation();
                    continue;
                }

                if (TreeDistS(T_work, T_end, child_k, empty_S)) {
                    debugPrint("TreeDistS: Found solution with rotation");
                    rollbackRotation();
                    return commitResult(true);
                }
                rollbackRotation();
            } catch (const SearchAbortException&) {
                throw;
            } catch (...) {
                rollbackRotation();
                continue;
            }
        }
        return commitResult(false);
    } else {
        //                       (a) include neither, (b) include eᵢ only, (c) include eᵢ′ only
        //                       (skip choices that pick a non-edge or conflict with independence)."
        if (g_profile.enabled) {
            g_profile.s_branch_calls++;
        }
        debugPrint("TreeDistS: Implementing complete S branching");

        dedupeSPairs(S_work);
        const int current_conflicts = cachedConflictCount(T_work, T_end);
        const bool must_drop_conflict_now = (k_work == current_conflicts);
        std::unordered_map<std::uint64_t, int> edge_gain_hint;
        std::unordered_map<std::uint64_t, int> edge_free_hint;
        std::unordered_map<std::uint64_t, int> edge_lb_hint;
        const bool use_scored_pair_order =
            pairHeuristicOrderingEnabled() && (must_drop_conflict_now || S_work.size() <= 24);
        if (use_scored_pair_order) {
            scoreAndOrderSPairs(T_work, T_end, current_conflicts, S_work,
                                edge_gain_hint, edge_free_hint, edge_lb_hint);
        }
        auto edgeGainHint = [&](const std::pair<int, int>& e) {
            auto it = edge_gain_hint.find(edgeScoreKey(e));
            if (it == edge_gain_hint.end()) return std::numeric_limits<int>::min() / 2;
            return it->second;
        };
        auto edgeFreeHint = [&](const std::pair<int, int>& e) {
            auto it = edge_free_hint.find(edgeScoreKey(e));
            if (it == edge_free_hint.end()) return 0;
            return it->second;
        };
        auto edgeLbHint = [&](const std::pair<int, int>& e) {
            auto it = edge_lb_hint.find(edgeScoreKey(e));
            if (it == edge_lb_hint.end()) return std::numeric_limits<int>::max() / 2;
            return it->second;
        };
        auto edgeExistsNow = [&](const std::pair<int, int>& e) -> bool {
            auto oriented = orientEdge(T_work, e);
            return oriented.first >= 0 && oriented.second >= 0;
        };

        auto runBranchingWithPairs =
            [&](const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& active_pairs) -> bool {
            std::vector<std::pair<int, int>> chosen;
            std::unordered_set<int> used_nodes;
            bool use_mask_memo = (T_work.original_nodes.size() <= 63);
            std::unordered_map<int, int> node_to_bit;
            if (use_mask_memo) {
                int bit = 0;
                for (int node : T_work.original_nodes) {
                    node_to_bit[node] = bit++;
                }
            }
            auto edgeMask = [&](const std::pair<int, int>& e, uint64_t& mask_out) -> bool {
                if (!use_mask_memo) return false;
                auto it1 = node_to_bit.find(e.first);
                auto it2 = node_to_bit.find(e.second);
                if (it1 == node_to_bit.end() || it2 == node_to_bit.end()) return false;
                if (it1->second >= 63 || it2->second >= 63) return false;
                mask_out = (1ULL << it1->second) | (1ULL << it2->second);
                return true;
            };
            uint64_t used_mask = 0;
            struct BranchStateKey {
                uint32_t idx = 0;
                uint64_t mask = 0;
                bool operator==(const BranchStateKey& o) const {
                    return idx == o.idx && mask == o.mask;
                }
            };
            struct BranchStateHash {
                size_t operator()(const BranchStateKey& k) const {
                    uint64_t h = (static_cast<uint64_t>(k.idx) << 32) ^ k.mask;
                    h ^= (h >> 33);
                    h *= 0xff51afd7ed558ccdULL;
                    h ^= (h >> 33);
                    return static_cast<size_t>(h);
                }
            };
            std::unordered_set<BranchStateKey, BranchStateHash> failed_states;

            auto canInclude = [&](const std::pair<int, int>& e) -> bool {
                auto oriented = orientEdge(T_work, e);
                if (oriented.first < 0 || oriented.second < 0) {
                    return false;
                }
                if (use_mask_memo) {
                    uint64_t m = 0;
                    if (!edgeMask(e, m)) return false;
                    return (used_mask & m) == 0;
                }
                return used_nodes.count(e.first) == 0 && used_nodes.count(e.second) == 0;
            };

            std::function<bool(size_t)> branchPairs = [&](size_t idx) -> bool {
                BranchStateKey state{static_cast<uint32_t>(idx), used_mask};
                if (use_mask_memo && failed_states.find(state) != failed_states.end()) {
                    return false;
                }
                if ((int)chosen.size() > k_work) {
                    if (use_mask_memo) failed_states.insert(state);
                    return false;
                }
                if (idx == active_pairs.size()) {
                    if (chosen.empty()) {
                        if (use_mask_memo) failed_states.insert(state);
                        return false;
                    }
                    const Key128 i_bounds_key = makeKeyIBase(T_work, T_end, chosen);
                    bool i_bounds_value = false;
                    if (tryBoundsPrune(g_tdi_bounds, i_bounds_key, k_work, i_bounds_value)) {
                        if (!i_bounds_value && use_mask_memo) failed_states.insert(state);
                        return i_bounds_value;
                    }
                    int i_lb = combinedLowerBoundI(T_work, T_end, chosen, current_conflicts);
                    if (i_lb > k_work) {
                        if (use_mask_memo) failed_states.insert(state);
                        return false;
                    }
                    bool i_ok = TreeDistI(T_work, T_end, k_work, chosen);
                    if (!i_ok && use_mask_memo) failed_states.insert(state);
                    return i_ok;
                }

                const auto& pair = active_pairs[idx];
                const auto& e1 = pair.first;
                const auto& e2 = pair.second;
                if (!use_scored_pair_order) {
                    // Baseline Li–Xia branch order: include neither, include e1, include e2
                    if (branchPairs(idx + 1)) return true;

                    if (canInclude(e1)) {
                        uint64_t m1 = 0;
                        bool has_m1 = edgeMask(e1, m1);
                        chosen.push_back(e1);
                        used_nodes.insert(e1.first);
                        used_nodes.insert(e1.second);
                        if (has_m1) used_mask |= m1;
                        if (branchPairs(idx + 1)) return true;
                        if (has_m1) used_mask &= ~m1;
                        used_nodes.erase(e1.first);
                        used_nodes.erase(e1.second);
                        chosen.pop_back();
                    }

                    if (!(e1 == e2) && canInclude(e2)) {
                        uint64_t m2 = 0;
                        bool has_m2 = edgeMask(e2, m2);
                        chosen.push_back(e2);
                        used_nodes.insert(e2.first);
                        used_nodes.insert(e2.second);
                        if (has_m2) used_mask |= m2;
                        if (branchPairs(idx + 1)) return true;
                        if (has_m2) used_mask &= ~m2;
                        used_nodes.erase(e2.first);
                        used_nodes.erase(e2.second);
                        chosen.pop_back();
                    }

                    if (use_mask_memo) failed_states.insert(state);
                    return false;
                }

                auto tryIncludeEdge = [&](const std::pair<int, int>& e) -> bool {
                    if (!canInclude(e)) return false;
                    uint64_t m = 0;
                    bool has_m = edgeMask(e, m);
                    chosen.push_back(e);
                    used_nodes.insert(e.first);
                    used_nodes.insert(e.second);
                    if (has_m) used_mask |= m;
                    if (branchPairs(idx + 1)) return true;
                    if (has_m) used_mask &= ~m;
                    used_nodes.erase(e.first);
                    used_nodes.erase(e.second);
                    chosen.pop_back();
                    return false;
                };

                auto tryIncludeByPromise = [&]() -> bool {
                    std::pair<int, int> first = e1;
                    std::pair<int, int> second = e2;
                    int g1 = edgeGainHint(first);
                    int g2 = edgeGainHint(second);
                    int f1 = edgeFreeHint(first);
                    int f2 = edgeFreeHint(second);
                    int lb1 = edgeLbHint(first);
                    int lb2 = edgeLbHint(second);
                    if (lb2 < lb1 ||
                        (lb2 == lb1 && (g2 > g1 || (g2 == g1 && f2 > f1)))) {
                        std::swap(first, second);
                    }
                    if (tryIncludeEdge(first)) return true;
                    if (!(first == second) && tryIncludeEdge(second)) return true;
                    return false;
                };

                auto tryOptionsByFScore = [&]() -> bool {
                    struct Option {
                        int type = 0;  // 0: neither, 1: e1, 2: e2
                        int fscore = std::numeric_limits<int>::max() / 2;
                    };
                    std::vector<Option> options;
                    options.reserve(3);

                    int neither_lb = combinedLowerBoundI(T_work, T_end, chosen, current_conflicts);
                    options.push_back({0, neither_lb});

                    if (canInclude(e1)) {
                        std::vector<std::pair<int, int>> tmp = chosen;
                        tmp.push_back(e1);
                        int lb = combinedLowerBoundI(T_work, T_end, tmp, current_conflicts);
                        options.push_back({1, 1 + lb});
                    }
                    if (!(e1 == e2) && canInclude(e2)) {
                        std::vector<std::pair<int, int>> tmp = chosen;
                        tmp.push_back(e2);
                        int lb = combinedLowerBoundI(T_work, T_end, tmp, current_conflicts);
                        options.push_back({2, 1 + lb});
                    }

                    std::sort(options.begin(), options.end(), [](const Option& a, const Option& b) {
                        if (a.fscore != b.fscore) return a.fscore < b.fscore;
                        return a.type < b.type;
                    });

                    for (const auto& option : options) {
                        if (option.type == 0) {
                            if (branchPairs(idx + 1)) return true;
                        } else if (option.type == 1) {
                            if (tryIncludeEdge(e1)) return true;
                        } else if (option.type == 2) {
                            if (tryIncludeEdge(e2)) return true;
                        }
                    }
                    return false;
                };

                int best_gain = std::max(edgeGainHint(e1), edgeGainHint(e2));
                if (use_scored_pair_order) {
                    if (tryOptionsByFScore()) return true;
                } else {
                    bool include_first = must_drop_conflict_now || best_gain > 0;
                    if (include_first) {
                        if (tryIncludeByPromise()) return true;
                        if (branchPairs(idx + 1)) return true;
                    } else {
                        if (branchPairs(idx + 1)) return true;
                        if (tryIncludeByPromise()) return true;
                    }
                }

                if (use_mask_memo) failed_states.insert(state);
                return false;
            };

            return branchPairs(0);
        };

        std::vector<size_t> prefix_sizes;
        const size_t full_size = S_work.size();
        const size_t progressive_threshold =
            static_cast<size_t>(std::max(progressiveSPrefixStart(), 24));
        const bool use_progressive_prefix =
            progressiveSPrefixEnabled() && use_scored_pair_order && full_size > progressive_threshold;
        if (use_progressive_prefix && full_size > 0) {
            size_t cap = std::min(full_size, static_cast<size_t>(progressiveSPrefixStart()));
            prefix_sizes.push_back(cap);
            while (prefix_sizes.back() < full_size) {
                size_t next = std::min(full_size, prefix_sizes.back() * 2);
                if (next == prefix_sizes.back()) {
                    next = full_size;
                }
                prefix_sizes.push_back(next);
            }
        } else {
            prefix_sizes.push_back(full_size);
        }
        if (prefix_sizes.empty() || prefix_sizes.back() != full_size) {
            prefix_sizes.push_back(full_size);
        }

        for (size_t prefix : prefix_sizes) {
            std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> active_pairs;
            active_pairs.reserve(prefix);
            for (auto it = S_work.begin(); it != S_work.begin() + prefix; ++it) {
                const bool e1_live = edgeExistsNow(it->first);
                const bool e2_live = (!(it->first == it->second) && edgeExistsNow(it->second));
                if (!e1_live && !e2_live) {
                    if (g_profile.enabled) {
                        g_profile.s_dead_pair_skips++;
                    }
                    continue;
                }
                active_pairs.push_back(*it);
            }
            if (active_pairs.empty()) {
                if (g_profile.enabled) {
                    g_profile.s_dead_prefixes++;
                }
                if (prefix == full_size) {
                    if (g_profile.enabled) {
                        g_profile.s_full_dead_prefix_failures++;
                    }
                    return commitResult(false);
                }
                continue;
            }
            if (runBranchingWithPairs(active_pairs)) {
                debugPrint("TreeDistS: Found solution in S branching");
                return commitResult(true);
            }
        }
    }

    debugPrint("TreeDistS: No solution found");
    return commitResult(false);
}

// Definition: Find the smallest k such that FlipDistTree returns true
// Parameters: T_init/T_final: trees; max_k: search cap
// Returns: minimum k, or -1 if not found up to max_k
// Errors: returns -1 on invalid inputs or failures during search
int FlipDistMinK(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int max_k) {
    if (max_k < 0) return -1;
    const size_t n_for_reserve = T_init.original_nodes.size();
    g_active_problem_size = n_for_reserve;
    if (n_for_reserve >= 23) {
        const size_t state_reserve = n_for_reserve >= 25 ? 1200000ULL : (n_for_reserve == 24 ? 800000ULL : 650000ULL);
        const size_t side_reserve = state_reserve / 2;
        g_memo_tds.max_load_factor(0.5f);
        g_tds_bounds.max_load_factor(0.5f);
        g_conflict_cache.max_load_factor(0.5f);
        g_free_edge_cache.max_load_factor(0.5f);
        g_empty_s_child_cache.max_load_factor(0.5f);
        g_free_edge_partition_side_budget_cache.max_load_factor(0.5f);
        g_pair_incumbent_upper.max_load_factor(0.5f);
        g_memo_tds.reserve(state_reserve);
        g_tds_bounds.reserve(state_reserve);
        g_conflict_cache.reserve(state_reserve);
        g_free_edge_cache.reserve(state_reserve);
        g_empty_s_child_cache.reserve(side_reserve);
        g_free_edge_partition_side_budget_cache.reserve(side_reserve);
        g_pair_incumbent_upper.reserve(side_reserve);
    }
    const Key128 pair_key = makeKeyPair(T_init, T_final);
    KBounds &bounds = g_kbounds[pair_key];
    if (bounds.min_success >= 0 && bounds.max_fail >= bounds.min_success) {
        bounds.max_fail = bounds.min_success - 1;
    }

    if (bounds.min_success >= 0 && bounds.max_fail >= bounds.min_success - 1 &&
        bounds.min_success <= max_k) {
        return bounds.min_success;
    }

    int lo = std::max(0, bounds.max_fail + 1);
    // Exact lower bound: every conflicting edge must be removed at least once
    lo = std::max(lo, countConflictEdges(T_init, T_final));
    if (lo > max_k) {
        return -1;
    }

    int hi = -1;
    if (bounds.min_success >= 0 && bounds.min_success <= max_k) {
        hi = bounds.min_success;
    }

    // Find a feasible upper bound progressively instead of probing max_k directly
    if (hi < 0) {
        int probe = lo;
        const int linear_prefix_end = std::min(max_k, lo + minKLinearPrefixProbes(T_init, T_final));
        const int linear_abort_ms = minKLinearAbortMs(T_init);
        while (true) {
            const auto probe_t0 = std::chrono::steady_clock::now();
            const bool probe_ok = FlipDistTree(T_init, T_final, probe);
            const auto probe_t1 = std::chrono::steady_clock::now();
            const double probe_ms =
                std::chrono::duration<double, std::milli>(probe_t1 - probe_t0).count();
            if (probe_ok) {
                hi = probe;
                if (bounds.min_success < 0 || probe < bounds.min_success) {
                    bounds.min_success = probe;
                }
                break;
            }
            if (probe > bounds.max_fail) bounds.max_fail = probe;
            if (probe >= max_k) {
                return -1;
            }

            const int probe_step = minKProbeStep(T_init);
            const bool continue_linear =
                probe < linear_prefix_end &&
                !(linear_abort_ms > 0 && probe_ms >= static_cast<double>(linear_abort_ms));
            long long next_probe = continue_linear
                                       ? static_cast<long long>(probe) + 1LL
                                       : static_cast<long long>(probe) + static_cast<long long>(probe_step);
            if (next_probe > max_k) next_probe = max_k;
            if (next_probe <= probe) next_probe = static_cast<long long>(probe) + 1LL;
            probe = static_cast<int>(next_probe);
        }
    }

    if (lo > hi) {
        if (bounds.min_success >= 0 && bounds.min_success <= max_k) {
            return bounds.min_success;
        }
        return -1;
    }

    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (FlipDistTree(T_init, T_final, mid)) {
            hi = mid;
            if (bounds.min_success < 0 || mid < bounds.min_success) {
                bounds.min_success = mid;
            }
        } else {
            if (mid > bounds.max_fail) bounds.max_fail = mid;
            lo = mid + 1;
        }
    }

    if (bounds.min_success < 0 || lo < bounds.min_success) {
        bounds.min_success = lo;
    }
    if (lo - 1 > bounds.max_fail) bounds.max_fail = lo - 1;
    return lo;
}

void resetAlgorithmCaches() {
    g_conflict_cache.clear();
    g_target_edge_set_cache.clear();
    g_free_edge_cache.clear();
    g_heuristic_bundle_cache.clear();
    g_small_exact_dist_cache.clear();
    g_local_nopath_lb_cache.clear();
    g_tds_inflight_best_budget.clear();
    g_tdi_inflight_best_budget.clear();
    g_pair_incumbent_upper.clear();
    g_empty_s_core_cache.clear();
    g_empty_s_core_exact_cache.clear();
    g_empty_s_motif_summary_cache.clear();
    g_empty_s_motif_exact_cache.clear();
    g_branchy_core_summary_cache.clear();
    g_branchy_core_exact_cache.clear();
    g_pair_common_edge_cache.clear();
    g_empty_s_free_hint_gate_cache.clear();
    g_prefix_pair_request_memo.clear();
    g_empty_s_child_cache.clear();
    g_common_edge_decompose_cache.clear();
    g_common_edge_double_rotate_cache.clear();
    g_free_edge_partition_structure_cache.clear();
    g_target_partition_range_cache.clear();
    g_free_edge_partition_side_cache.clear();
    g_free_edge_partition_split_cache.clear();
    g_free_edge_partition_split_signature_cache.clear();
    g_free_edge_partition_split_budget_cache.clear();
    g_free_edge_partition_side_budget_cache.clear();
    g_free_edge_partition_side_budget_range_cache.clear();
    g_plateau_tree_dumped_states.clear();
    g_plateau_tree_dump_count = 0;
}
