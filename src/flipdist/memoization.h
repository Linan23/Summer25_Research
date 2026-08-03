#pragma once

#include "bf_bst.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef FLIPDIST_REUSE_PROFILE_COUNTERS
#define FLIPDIST_REUSE_PROFILE_COUNTERS 0
#endif

inline constexpr bool DEBUG = false;

struct ProfileStats {
    bool enabled = false;
    long long calls_flipdist = 0;
    long long calls_tdi = 0;
    long long calls_tds = 0;
    long long tdi_bounds_prunes = 0;
    long long tdi_conflict_lb_prunes = 0;
    long long tdi_remaining_budget_prunes = 0;
    long long tdi_exact_budget_success = 0;
    long long tdi_exact_budget_fail = 0;
    long long tdi_tree_equal_success = 0;
    long long tdi_post_i_empty_s = 0;
    long long tdi_post_i_nonempty_s = 0;
    long long tdi_post_i_lb_prunes = 0;
    long long tdi_prefix_lb_checks = 0;
    long long tdi_prefix_lb_prunes = 0;
    long long tdi_prefix_pairgen_calls = 0;
    long long tdi_prefix_pairgen_pairs_scanned = 0;
    long long tdi_prefix_pairgen_pairs_inserted = 0;
    long long tdi_prefix_request_total = 0;
    long long tdi_prefix_request_unique = 0;
    long long tdi_prefix_request_repeated = 0;
    long long tdi_prefix_emitted_pairs_0 = 0;
    long long tdi_prefix_emitted_pairs_1 = 0;
    long long tdi_prefix_emitted_pairs_2 = 0;
    long long tdi_prefix_conflict_calls = 0;
    long long tdi_prefix_bound_calls = 0;
    long long tdi_tds_success = 0;
    long long tdi_tds_fail = 0;
    long long tdi_rotation_skips = 0;
    long long memo_hits_tds = 0;
    long long bounds_hits_tds = 0;
    long long tds_unique_states = 0;
    long long tds_repeated_states = 0;
    long long calls_partition = 0;
    long long partition_unique_structures = 0;
    long long partition_repeated_structures = 0;
    long long partition_unique_side_states = 0;
    long long partition_repeated_side_states = 0;
    long long partition_unique_split_signatures = 0;
    long long partition_repeated_split_signatures = 0;
    long long partition_precheck_calls = 0;
    long long partition_precheck_matches = 0;
    long long partition_precheck_mismatches = 0;
    long long partition_budget_iterations = 0;
    long long tds_partition_side_budget_cache_hits = 0;
    long long tds_partition_side_budget_cache_misses = 0;
    long long tds_partition_split_cache_hits = 0;
    long long tds_partition_split_cache_misses = 0;
    long long partition_empty_side_shortcuts = 0;
    long long partition_side1_recursions = 0;
    long long partition_side2_recursions = 0;
    long long partition_side1_bound_hits = 0;
    long long partition_side2_bound_hits = 0;
    long long free_edge_hits = 0;
    long long free_edge_misses = 0;
    long long s_branch_calls = 0;
    long long s_empty_branch_calls = 0;
    long long s_empty_unique_states = 0;
    long long s_empty_repeated_states = 0;
    long long s_empty_no_candidate_states = 0;
    long long s_empty_progress_states = 0;
    long long s_empty_plateau_states = 0;
    long long s_empty_progress_candidates = 0;
    long long s_empty_plateau_candidates = 0;
    long long s_empty_plateau_buckets = 0;
    long long s_empty_artic_arm_hits = 0;
    long long s_empty_duplicate_child_states = 0;
    long long s_empty_must_drop_rejects = 0;
    long long s_empty_incumbent_rejects = 0;
    long long s_empty_core_lookups = 0;
    long long s_empty_core_cache_hits = 0;
    long long s_empty_core_exact_hits = 0;
    long long s_empty_core_common_splits = 0;
    long long s_empty_core_free_reductions = 0;
    long long s_empty_motif_common_peels = 0;
    long long s_empty_motif_free_peels = 0;
    long long s_empty_motif_exact_fallbacks = 0;
    long long s_empty_motif_exact_cache_hits = 0;
    long long s_empty_branchy_reduction_hits = 0;
    long long s_empty_branchy_exact_bounds = 0;
    long long s_empty_branchy_exact_decisions = 0;
    long long s_empty_branchy_exact_cache_hits = 0;
    long long s_dead_pair_skips = 0;
    long long s_dead_prefixes = 0;
    long long s_full_dead_prefix_failures = 0;
    long long independent_subsets_initial = 0;
    long long independent_subsets_s = 0;
    long long max_s_size = 0;
    long long max_i_size = 0;
    long long heuristic_lb_calls = 0;
    long long heuristic_lb_cache_hits = 0;
    long long conflict_cache_hits = 0;
    long long free_edge_cache_hits = 0;
    long long dominance_prunes = 0;
    long long small_exact_lookups = 0;
    long long small_exact_hits = 0;
    long long local_nopath_lookups = 0;
    long long local_nopath_hits = 0;
    double time_flipdist_ms = 0.0;
    double time_tdi_ms = 0.0;
    double time_tdi_prefix_pairgen_ms = 0.0;
    double time_tdi_prefix_pairgen_helper_ms = 0.0;
    double time_tdi_prefix_pairgen_build_ms = 0.0;
    double time_tdi_prefix_pairgen_build_edges_ms = 0.0;
    double time_tdi_prefix_pairgen_build_ranges_ms = 0.0;
    double time_tdi_prefix_pairgen_build_range_lookup_ms = 0.0;
    double time_tdi_prefix_pairgen_build_insert_ms = 0.0;
    double time_tdi_prefix_pairgen_build_boundary_ms = 0.0;
    double time_tdi_prefix_pairgen_pick_ms = 0.0;
    double time_tdi_prefix_pairgen_insert_ms = 0.0;
    double time_tdi_prefix_conflict_ms = 0.0;
    double time_tdi_prefix_bound_ms = 0.0;
    double time_tds_ms = 0.0;
    double time_tds_free_filter_ms = 0.0;
    double time_tds_free_rotate_ms = 0.0;
    double time_tds_free_partner_ms = 0.0;
    double time_partition_ms = 0.0;
    double time_partition_precheck_ms = 0.0;
    double time_partition_build_ms = 0.0;
    double time_partition_side_stats_ms = 0.0;
    double time_partition_count_edges_ms = 0.0;
    double time_partition_conflicts_ms = 0.0;
    double time_partition_bounds_ms = 0.0;
    double time_partition_split_s_ms = 0.0;
    double time_partition_budget_loop_ms = 0.0;
    std::vector<int> tdi_prefix_request_vcount_samples;
    std::vector<int> tdi_prefix_request_bucket_l_samples;
    std::vector<int> tdi_prefix_request_bucket_r_samples;
    std::chrono::steady_clock::time_point start_time;
    int abort_ms = -1;
};

extern ProfileStats g_profile;

struct Key128 {
    std::uint64_t hi = 0;
    std::uint64_t lo = 0;
    // Purpose: Compares two keys or state signatures for exact equality.
    // Inputs: o: the value to compare with this instance.
    // Output: True when all identity fields match; otherwise false.
    bool operator==(const Key128 &o) const {
        return hi == o.hi && lo == o.lo;
    }
};

struct Key128Hash {
    // Purpose: Hashes a key or state signature for unordered containers.
    // Inputs: The key or signature supplied to the call operator.
    // Output: A hash code consistent with the corresponding equality operator.
    std::size_t operator()(const Key128 &k) const;
};

struct PlateauStateRecord {
    Key128 state_key{};
    int internal_edges = 0;
    int conflicts = 0;
    int reduced_core_internal_edges = 0;
    int min_k = -1;
    int max_k = -1;
    int legal_children = -1;
    int plateau_buckets = -1;
    int start_branching_nodes = 0;
    int target_branching_nodes = 0;
    int start_max_branch_subtree_edges = 0;
    int target_max_branch_subtree_edges = 0;
    int chain_arm_count = 0;
    int broom_arm_count = 0;
    int branchy_reduction_hits = 0;
    int branchy_forced_cost_max = 0;
    int branchy_residual_edges_min = -1;
    int branchy_exact_bound_hits = 0;
    int branchy_exact_decision_hits = 0;
    int recurrence_count = 0;
    int success_count = 0;
    int fail_count = 0;
    int timeout_count = 0;
    int free_edge_hits = 0;
    int free_edge_misses = 0;
    std::string motif = "unknown";
    std::string core_signature;
    std::string branchy_residual_motif;
    std::string start_degree_histogram;
    std::string target_degree_histogram;
};

struct TdiPostIStateRecord {
    Key128 state_key{};
    int conflicts = 0;
    int s_size = 0;
    int min_remaining_k = -1;
    int max_remaining_k = -1;
    int min_lower_bound = -1;
    int max_lower_bound = -1;
    int start_branching_nodes = 0;
    int target_branching_nodes = 0;
    int recurrence_count = 0;
    std::string start_degree_histogram;
    std::string target_degree_histogram;
    std::string start_tree;
    std::string target_tree;
};

struct KBounds {
    int max_fail = -1;
    int min_success = -1;
};

extern std::unordered_map<Key128, bool, Key128Hash> g_memo_flipdist;
extern std::unordered_map<Key128, bool, Key128Hash> g_memo_tdi;
extern std::unordered_map<Key128, bool, Key128Hash> g_memo_tds;
extern std::unordered_map<Key128, KBounds, Key128Hash> g_kbounds;
extern std::unordered_map<Key128, KBounds, Key128Hash> g_tds_bounds;
extern std::unordered_map<Key128, KBounds, Key128Hash> g_tdi_bounds;

// Purpose: Reports whether the plateau profile feature is enabled.
// Inputs: None.
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool plateauProfileEnabled();
// Purpose: Resets plateau profile.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void resetPlateauProfile();
// Purpose: Starts plateau profile run.
// Inputs: case_type (const std::string &); n (int); seed (int); direction (const std::string &); clear_output_file (bool).
// Output: No return value; updates referenced outputs or shared state.
void beginPlateauProfileRun(const std::string &case_type, int n, int seed,
                            const std::string &direction, bool clear_output_file);
// Purpose: Records plateau state sample.
// Inputs: state_key (const Key128 &); internal_edges (int); reduced_core_internal_edges (int); k (int); conflicts (int); legal_children (int); plateau_buckets (int); free_edge_hit (bool); start_branching_nodes (int); target_branching_nodes (int); start_max_branch_subtree_edges (int); target_max_branch_subtree_edges (int); chain_arm_count (int); broom_arm_count (int); start_degree_histogram (const std::string &); target_degree_histogram (const std::string &); motif (const std::string &); core_signature (const std::string &); branchy_reduction_triggered (bool); branchy_forced_cost (int); branchy_residual_edges (int); branchy_residual_motif (const std::string &); branchy_exact_bound_used (bool); branchy_exact_decision_used (bool).
// Output: No return value; updates referenced outputs or shared state.
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
                            bool branchy_exact_bound_used, bool branchy_exact_decision_used);
// Purpose: Records plateau state outcome.
// Inputs: state_key (const Key128 &); outcome (const std::string &).
// Output: No return value; updates referenced outputs or shared state.
void notePlateauStateOutcome(const Key128 &state_key, const std::string &outcome);
// Purpose: Finalizes plateau profile run.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void finalizePlateauProfileRun();
// Purpose: Emits plateau profile jsonl.
// Inputs: os (std::ostream &); case_type (const std::string &); n (int); seed (int); direction (const std::string &).
// Output: No return value; updates referenced outputs or shared state.
void emitPlateauProfileJsonl(std::ostream &os, const std::string &case_type, int n, int seed,
                             const std::string &direction);
// Purpose: Reports whether the tdi post i profile feature is enabled.
// Inputs: None.
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool tdiPostIProfileEnabled();
// Purpose: Resets tdi post i profile.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void resetTdiPostIProfile();
// Purpose: Starts tdi post i profile run.
// Inputs: case_type (const std::string &); n (int); seed (int); direction (const std::string &); clear_output_file (bool).
// Output: No return value; updates referenced outputs or shared state.
void beginTdiPostIProfileRun(const std::string &case_type, int n, int seed,
                             const std::string &direction, bool clear_output_file);
// Purpose: Records tdi post i prune sample.
// Inputs: state_key (const Key128 &); remaining_k (int); lower_bound (int); conflicts (int); s_size (int); start_branching_nodes (int); target_branching_nodes (int); start_degree_histogram (const std::string &); target_degree_histogram (const std::string &); start_tree (const std::string &); target_tree (const std::string &).
// Output: No return value; updates referenced outputs or shared state.
void noteTdiPostIPruneSample(const Key128 &state_key, int remaining_k, int lower_bound,
                             int conflicts, int s_size, int start_branching_nodes,
                             int target_branching_nodes,
                             const std::string &start_degree_histogram,
                             const std::string &target_degree_histogram,
                             const std::string &start_tree,
                             const std::string &target_tree);
// Purpose: Emits tdi post i profile jsonl.
// Inputs: os (std::ostream &); case_type (const std::string &); n (int); seed (int); direction (const std::string &).
// Output: No return value; updates referenced outputs or shared state.
void emitTdiPostIProfileJsonl(std::ostream &os, const std::string &case_type, int n, int seed,
                              const std::string &direction);
// Purpose: Finalizes tdi post i profile run.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void finalizeTdiPostIProfileRun();

// Purpose: Removes duplicates from s pairs.
// Inputs: pairs (std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &).
// Output: No return value; updates referenced outputs or shared state.
void dedupeSPairs(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &pairs);

// Purpose: Builds key pair.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &).
// Output: The resulting Key128 value.
Key128 makeKeyPair(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);
// Purpose: Builds key flip.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); k (int).
// Output: The resulting Key128 value.
Key128 makeKeyFlip(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k);
// Purpose: Builds key i.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); k (int); I (const std::vector<std::pair<int, int>> &).
// Output: The resulting Key128 value.
Key128 makeKeyI(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k,
                const std::vector<std::pair<int, int>> &I);
// Purpose: Builds key s.
// Inputs: T_init (const VectorRangeTreeMap &); T_end (const VectorRangeTreeMap &); k (int); S (const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &).
// Output: The resulting Key128 value.
Key128 makeKeyS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
                const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);
// Purpose: Builds key s base.
// Inputs: T_init (const VectorRangeTreeMap &); T_end (const VectorRangeTreeMap &); S (const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &).
// Output: The resulting Key128 value.
Key128 makeKeySBase(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end,
                    const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);
// Purpose: Builds key i base.
// Inputs: T_init (const VectorRangeTreeMap &); T_final (const VectorRangeTreeMap &); I (const std::vector<std::pair<int, int>> &).
// Output: The resulting Key128 value.
Key128 makeKeyIBase(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final,
                    const std::vector<std::pair<int, int>> &I);

struct EmptySKeys {
    Key128 pair_key{};
    Key128 bounds_key{};
    Key128 exact_key{};
    Key128 exact_base_key{};
};

struct EmptySBaseKeys {
    Key128 pair_key{};
    Key128 bounds_key{};
};

// Purpose: Builds empty s base keys.
// Inputs: T_init (const VectorRangeTreeMap &); T_end (const VectorRangeTreeMap &).
// Output: The resulting EmptySBaseKeys value.
EmptySBaseKeys makeEmptySBaseKeys(const VectorRangeTreeMap &T_init,
                                  const VectorRangeTreeMap &T_end);
// Purpose: Builds empty s keys.
// Inputs: T_init (const VectorRangeTreeMap &); T_end (const VectorRangeTreeMap &); k (int).
// Output: The resulting EmptySKeys value.
EmptySKeys makeEmptySKeys(const VectorRangeTreeMap &T_init,
                          const VectorRangeTreeMap &T_end,
                          int k);
// Purpose: Builds empty s exact key from base.
// Inputs: exact_base_key (Key128); k (int).
// Output: The resulting Key128 value.
Key128 makeEmptySExactKeyFromBase(Key128 exact_base_key, int k);

// Purpose: Attempts bounds prune.
// Inputs: bounds_map (const std::unordered_map<Key128, KBounds, Key128Hash> &); key (const Key128 &); k (int); value_out (bool &).
// Output: True when the stated condition holds or operation succeeds; otherwise false.
inline bool tryBoundsPrune(const std::unordered_map<Key128, KBounds, Key128Hash> &bounds_map,
                           const Key128 &key, int k, bool &value_out) {
    auto it = bounds_map.find(key);
    if (it == bounds_map.end()) return false;
    const auto &b = it->second;
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

// Purpose: Updates bounds.
// Inputs: bounds_map (std::unordered_map<Key128, KBounds, Key128Hash> &); key (const Key128 &); k (int); value (bool).
// Output: No return value; updates referenced outputs or shared state.
inline void updateBounds(std::unordered_map<Key128, KBounds, Key128Hash> &bounds_map,
                         const Key128 &key, int k, bool value) {
    auto &b = bounds_map[key];
    if (value) {
        if (b.min_success < 0 || k < b.min_success) {
            b.min_success = k;
        }
        const int forced_fail = k - 1;
        if (forced_fail < b.max_fail) {
            b.max_fail = forced_fail;
        }
        if (b.max_fail >= b.min_success) {
            b.max_fail = b.min_success - 1;
        }
    } else if (k > b.max_fail) {
        b.max_fail = k;
    }
}

// Purpose: Derives the first unfailed budget from learned bounds.
// Inputs: bounds_map (const std::unordered_map<Key128, KBounds, Key128Hash> &); key (const Key128 &).
// Output: The computed value or status code.
inline int requiredBudgetFromBounds(const std::unordered_map<Key128, KBounds, Key128Hash> &bounds_map,
                                    const Key128 &key) {
    auto it = bounds_map.find(key);
    if (it == bounds_map.end()) return 0;
    return std::max(0, it->second.max_fail + 1);
}

// Purpose: Resets memo.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void resetMemo();
// Purpose: Initializes profiling state from the active environment configuration.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void initProfile();
// Purpose: Reports abort requested.
// Inputs: None.
// Output: True when the stated condition holds or operation succeeds; otherwise false.
bool profileAbortRequested();
struct SearchAbortException {};
// Purpose: Throws when if profile abort.
// Inputs: None.
// Output: No return value; updates referenced outputs or shared state.
void throwIfProfileAbort();

struct ScopedTimer {
    double *acc = nullptr;
    std::chrono::steady_clock::time_point t0;
    // Purpose: Starts timing when an accumulator is provided.
    // Inputs: a: optional millisecond accumulator to update when the scope ends.
    // Output: A scope timer whose destructor records elapsed time.
    explicit ScopedTimer(double *a) : acc(a) {
        if (acc) {
            t0 = std::chrono::steady_clock::now();
        }
    }
    // Purpose: Adds this scope's elapsed milliseconds to the configured accumulator.
    // Inputs: None.
    // Output: No return value; updates the accumulator passed to the constructor.
    ~ScopedTimer() {
        if (acc) {
            auto t1 = std::chrono::steady_clock::now();
            *acc += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
    }
};

// Purpose: Emits a diagnostic message when debug logging is enabled.
// Inputs: msg (const std::string &).
// Output: No return value; updates referenced outputs or shared state.
void debugPrint(const std::string &msg);
