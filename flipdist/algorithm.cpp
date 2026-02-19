#include "algorithm.h"
#include "helpers.h"
#include "memoization.h"

#include <algorithm>
#include <functional>

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
    if (profileAbortRequested()) {
        return false;
    }
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    const Key128 memo_key = makeKeyFlip(T_init, T_final, k);
    const Key128 pair_key = makeKeyPair(T_init, T_final);
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
        return value;
    };

    // BASE CASE: Check if trees are already identical
    if (TreesEqual(T_init, T_final)) {
        debugPrint("Trees already equal, returning true");
        return memoReturn(true);
    }

    bool handled_common = false;
    if (tryCommonEdgeDecomposition(T_init, T_final, k, handled_common)) {
        debugPrint("Common edge decomposition succeeded");
        return memoReturn(true);
    }
    if (handled_common) {
        debugPrint("Common edge decomposition failed");
        return memoReturn(false);
    }

    auto [hasFree, freeEdge] = findFreeEdge(T_init, T_final);
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
    if (profileAbortRequested()) {
        return false;
    }
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    const Key128 memo_key = makeKeyI(T_init, T_final, k, I);
    auto memo_it = g_memo_tdi.find(memo_key);
    if (memo_it != g_memo_tdi.end()) {
        return memo_it->second;
    }
    const Key128 bounds_key = makeKeyIBase(T_init, T_final, I);
    bool bounds_value = false;
    if (tryBoundsPrune(g_tdi_bounds, bounds_key, k, bounds_value)) {
        g_memo_tdi[memo_key] = bounds_value;
        return bounds_value;
    }
    auto memoReturn = [&](bool value) {
        g_memo_tdi[memo_key] = value;
        updateBounds(g_tdi_bounds, bounds_key, k, value);
        return value;
    };

    // Exact lower bound: each conflicting edge must be removed by at least one rotation
    int conflict_lb = countConflictEdges(T_init, T_final);
    if (conflict_lb > k) {
        return memoReturn(false);
    }

    int remaining_budget = k - (int)I.size();  // k - |I| from pseudocode

    if (remaining_budget < 0) {  // This covers φ(T_init) > k − |I| case
        debugPrint("TreeDistI: Not enough budget for |I| rotations");
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
            } catch (...) {
                continue;
            }
        }

        if (TreesEqual(T_test, T_final)) {
            debugPrint("TreeDistI: Solved with exactly |I| rotations");
            return memoReturn(true);
        } else {
            debugPrint("TreeDistI: No budget left and not solved");
            return memoReturn(false);
        }
    }

    debugPrint("TreeDistI: Proceeding with remaining_budget=" + std::to_string(remaining_budget));

    if (TreesEqual(T_init, T_final)) {
        debugPrint("TreeDistI: Trees already equal");
        return memoReturn(true);
    }

    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S;
    std::vector<std::pair<int,int>> newDiagonals;

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
            continue;
        }

        debugPrint("TreeDistI: Applied rotation, new edge ē connects " + std::to_string(u) + " and " + std::to_string(v));

        auto cr = T_bar.getRange(v);
        newDiagonals.emplace_back(cr.first, cr.second);
    }

    appendPartnerPairsFromDiagonals(T_bar, newDiagonals, S);

    int post_i_conflict_lb = countConflictEdges(T_bar, T_final);
    if (post_i_conflict_lb > k - (int)I.size()) {
        return memoReturn(false);
    }

    return memoReturn(TreeDistS(T_bar, T_final, k - (int)I.size(), S));
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
    if (profileAbortRequested()) {
        return false;
    }
    debugPrint("Entering TreeDistS with k=" + std::to_string(k) + ", |S|=" + std::to_string(S.size()));

    const Key128 memo_key = makeKeyS(T_init, T_end, k, S);
    auto memo_it = g_memo_tds.find(memo_key);
    if (memo_it != g_memo_tds.end()) {
        return memo_it->second;
    }
    const Key128 bounds_key = makeKeySBase(T_init, T_end, S);
    bool bounds_value = false;
    if (tryBoundsPrune(g_tds_bounds, bounds_key, k, bounds_value)) {
        g_memo_tds[memo_key] = bounds_value;
        return bounds_value;
    }
    auto memoReturn = [&](bool value) {
        g_memo_tds[memo_key] = value;
        updateBounds(g_tds_bounds, bounds_key, k, value);
        return value;
    };

    // Base case: trees already equal
    if (TreesEqual(T_init, T_end)) {
        debugPrint("TreeDistS: Trees already equal");
        return memoReturn(true);
    }

    if (k < 0) {
        debugPrint("TreeDistS: Negative budget");
        return memoReturn(false);
    }
    if (k == 0) {
        debugPrint("TreeDistS: No budget and trees differ");
        return memoReturn(false);
    }

    int conflict_lb = countConflictEdges(T_init, T_end);
    if (conflict_lb > k) {
        return memoReturn(false);
    }

    VectorRangeTreeMap T_work = safeCopyTree(T_init);
    int k_work = k;
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_work = S;
    dedupeSPairs(S_work);
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

    while (k_work > 0) {
        int work_conflict_lb = countConflictEdges(T_work, T_end);
        if (work_conflict_lb > k_work) {
            return commitResult(false);
        }

        const Key128 state_key = makeKeyS(T_work, T_end, k_work, S_work);
        auto state_it = g_memo_tds.find(state_key);
        if (state_it != g_memo_tds.end()) {
            return commitResult(state_it->second);
        }
        const Key128 state_bounds_key = makeKeySBase(T_work, T_end, S_work);
        bool state_bounds_value = false;
        if (tryBoundsPrune(g_tds_bounds, state_bounds_key, k_work, state_bounds_value)) {
            return commitResult(state_bounds_value);
        }
        unfolded_bounds_entries.push_back({state_key, state_bounds_key, k_work});

        auto [hasFree, freeEdge] = findFreeEdge(T_work, T_end);
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
            S_filtered.reserve(S_work.size() + 2);
            for (const auto& pair : S_work) {
                if (!(pair.first == freeEdge || pair.second == freeEdge)) {
                    S_filtered.push_back(pair);
                }
            }
            debugPrint("TreeDistS: Filtered S from " + std::to_string(S_work.size()) + " to " + std::to_string(S_filtered.size()) + " pairs");

            VectorRangeTreeMap T_bar = safeCopyTree(T_work);
            int u, v;

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

            if (TreesEqual(T_bar, T_end)) {
                debugPrint("TreeDistS: Solved with free edge rotation");
                return commitResult(true);
            }

            auto cr = T_bar.getRange(v);
            int L = cr.first;
            int R = cr.second;

            appendPartnerPairsFromDiagonals(T_bar, {{L, R}}, S_filtered);
            dedupeSPairs(S_filtered);

            auto parent_range = T_bar.getRange(u);
            auto child_range = T_bar.getRange(v);

            try {
                if (g_profile.enabled) {
                    g_profile.calls_partition++;
                }
                ScopedTimer _tp(g_profile.enabled ? &g_profile.time_partition_ms : nullptr);
                auto [T_bar1, T_bar2] = VectorRangeTreeMap::partitionAlongEdge(T_bar, parent_range, child_range);
                auto [T_end1, T_end2] = VectorRangeTreeMap::partitionAlongEdge(T_end, parent_range, child_range);

                if (T_bar1.original_nodes != T_end1.original_nodes ||
                    T_bar2.original_nodes != T_end2.original_nodes) {
                    debugPrint("TreeDistS: Partition mismatch, continuing free-edge contraction");
                    T_work = std::move(T_bar);
                    S_work = std::move(S_filtered);
                    --k_work;
                    continue;
                }

                auto [S1, S2] = partitionS(S_filtered, T_bar1, T_bar2);
                int n1 = countInternalEdges(T_bar1);
                int n2 = countInternalEdges(T_bar2);
                int conf1 = countConflictEdges(T_bar1, T_end1);
                int conf2 = countConflictEdges(T_bar2, T_end2);

                if (n1 == 0) {
                    if (conf2 > k_work - 1) return commitResult(false);
                    return commitResult(TreeDistS(T_bar2, T_end2, k_work - 1, S2));
                }

                if (n2 == 0) {
                    if (conf1 > k_work - 1) return commitResult(false);
                    return commitResult(TreeDistS(T_bar1, T_end1, k_work - 1, S1));
                }

                const Key128 s1_bounds_key = makeKeySBase(T_bar1, T_end1, S1);
                const Key128 s2_bounds_key = makeKeySBase(T_bar2, T_end2, S2);
                const int lb1 = std::max(requiredBudgetFromBounds(g_tds_bounds, s1_bounds_key), conf1);
                const int lb2 = std::max(requiredBudgetFromBounds(g_tds_bounds, s2_bounds_key), conf2);

                int min_k1 = std::max(0, lb1);
                int max_k1 = std::min(k_work - 1, k_work - 1 - lb2);
                if (min_k1 > max_k1) {
                    return commitResult(false);
                }

                for (int k1 = min_k1; k1 <= max_k1; k1++) {
                    int k2 = k_work - 1 - k1;

                    bool s1_known = false;
                    bool s1_ok = false;
                    s1_known = tryBoundsPrune(g_tds_bounds, s1_bounds_key, k1, s1_ok);
                    if (!s1_known) {
                        s1_ok = TreeDistS(T_bar1, T_end1, k1, S1);
                    }
                    if (!s1_ok) {
                        continue;
                    }

                    bool s2_known = false;
                    bool s2_ok = false;
                    s2_known = tryBoundsPrune(g_tds_bounds, s2_bounds_key, k2, s2_ok);
                    if (!s2_known) {
                        s2_ok = TreeDistS(T_bar2, T_end2, k2, S2);
                    }
                    if (s2_ok) {
                        return commitResult(true);
                    }
                }
                return commitResult(false);
            } catch (...) {
                debugPrint("TreeDistS: Partitioning failed, continuing free-edge contraction");
                T_work = std::move(T_bar);
                S_work = std::move(S_filtered);
                --k_work;
                continue;
            }
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
        const int current_conflicts = countConflictEdges(T_work, T_end);
        const bool must_drop_conflict_now = (k_work == current_conflicts);
        const int child_k = k_work - 1;
        const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> empty_S;
        std::unordered_set<Key128, Key128Hash> seen_child_states;
        auto edges = getInternalEdges(T_work);
        struct RotationCandidate {
            std::pair<int,int> edge;
            int next_conflicts = INT_MAX;
            int free_hint = 0;
            Key128 child_exact_key{};
            Key128 child_bounds_key{};
        };
        std::vector<RotationCandidate> candidates;
        candidates.reserve(edges.size());
        std::vector<RotationCandidate> plateau_candidates;
        plateau_candidates.reserve(edges.size());
        auto byPromise = [](const RotationCandidate& a, const RotationCandidate& b) {
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

        for (const auto& edge : edges) {
            try {
                VectorRangeTreeMap probe = safeCopyTree(T_work);
                int parent = edge.first;
                int child = edge.second;
                if (probe.getLeftChild(parent) == child) {
                    probe.rotateRight(parent);
                } else if (probe.getRightChild(parent) == child) {
                    probe.rotateLeft(parent);
                } else {
                    continue;
                }
                int next_conflicts = countConflictEdges(probe, T_end);
                if (next_conflicts > k_work - 1) {
                    continue;
                }
                Key128 child_exact_key = makeKeyS(probe, T_end, child_k, empty_S);
                if (!seen_child_states.insert(child_exact_key).second) {
                    continue;
                }
                auto child_memo_it = g_memo_tds.find(child_exact_key);
                if (child_memo_it != g_memo_tds.end()) {
                    if (child_memo_it->second) {
                        return commitResult(true);
                    }
                    continue;
                }
                Key128 child_bounds_key = makeKeySBase(probe, T_end, empty_S);
                bool known_tds = false;
                if (tryBoundsPrune(g_tds_bounds, child_bounds_key, child_k, known_tds)) {
                    if (known_tds) {
                        return commitResult(true);
                    }
                    continue;
                }
                // Exact lower-bound gate for this child state before recursing
                Key128 child_pair_key = makeKeyPair(probe, T_end);
                int child_lb = std::max(next_conflicts,
                                        std::max(requiredBudgetFromBounds(g_tds_bounds, child_bounds_key),
                                                 requiredBudgetFromBounds(g_kbounds, child_pair_key)));
                if (child_lb > child_k) {
                    continue;
                }
                auto [probe_has_free, _] = findFreeEdge(probe, T_end);
                RotationCandidate cand{edge, next_conflicts, probe_has_free ? 1 : 0,
                                       child_exact_key, child_bounds_key};
                if (must_drop_conflict_now && next_conflicts >= current_conflicts) {
                    continue;
                }
                if (next_conflicts < current_conflicts || cand.free_hint) {
                    candidates.push_back(cand);
                } else {
                    plateau_candidates.push_back(cand);
                }
            } catch (...) {
                continue;
            }
        }

        std::sort(candidates.begin(), candidates.end(), byPromise);
        std::sort(plateau_candidates.begin(), plateau_candidates.end(), byPromise);
        candidates.insert(candidates.end(), plateau_candidates.begin(), plateau_candidates.end());

        for (const auto& cand : candidates) {
            try {
                auto child_memo_it = g_memo_tds.find(cand.child_exact_key);
                if (child_memo_it != g_memo_tds.end()) {
                    if (child_memo_it->second) {
                        return commitResult(true);
                    }
                    continue;
                }
                bool known_tds = false;
                if (tryBoundsPrune(g_tds_bounds, cand.child_bounds_key, child_k, known_tds)) {
                    if (known_tds) {
                        return commitResult(true);
                    }
                    continue;
                }
                VectorRangeTreeMap T_test = safeCopyTree(T_work);
                int parent = cand.edge.first;
                int child = cand.edge.second;

                if (T_test.getLeftChild(parent) == child) {
                    T_test.rotateRight(parent);
                } else if (T_test.getRightChild(parent) == child) {
                    T_test.rotateLeft(parent);
                }

                if (TreeDistS(T_test, T_end, child_k, empty_S)) {
                    debugPrint("TreeDistS: Found solution with rotation");
                    return commitResult(true);
                }
            } catch (...) {
                continue;
            }
        }
    } else {
        //                       (a) include neither, (b) include eᵢ only, (c) include eᵢ′ only
        //                       (skip choices that pick a non-edge or conflict with independence)."
        if (g_profile.enabled) {
            g_profile.s_branch_calls++;
        }
        debugPrint("TreeDistS: Implementing complete S branching");

        std::vector<std::pair<int,int>> chosen;
        std::unordered_set<int> used_nodes;
        dedupeSPairs(S_work);

        bool use_mask_memo = (T_work.original_nodes.size() <= 63);
        std::unordered_map<int, int> node_to_bit;
        if (use_mask_memo) {
            int bit = 0;
            for (int node : T_work.original_nodes) {
                node_to_bit[node] = bit++;
            }
        }
        auto edgeMask = [&](const std::pair<int,int>& e, uint64_t& mask_out) -> bool {
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

        auto canInclude = [&](const std::pair<int,int>& e) -> bool {
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
            if (idx == S_work.size()) {
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
                bool i_ok = TreeDistI(T_work, T_end, k_work, chosen);
                if (!i_ok && use_mask_memo) failed_states.insert(state);
                return i_ok;
            }

            // Option 1: include neither
            if (branchPairs(idx + 1)) return true;

            const auto& pair = S_work[idx];
            const auto& e1 = pair.first;
            const auto& e2 = pair.second;

            // Option 2: include e1
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

            // Option 3: include e2 (avoid duplicate edge)
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
        };

        if (branchPairs(0)) {
            debugPrint("TreeDistS: Found solution in S branching");
            return commitResult(true);
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
    const Key128 pair_key = makeKeyPair(T_init, T_final);
    auto &bounds = g_kbounds[pair_key];
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
        while (true) {
            if (FlipDistTree(T_init, T_final, probe)) {
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

            long long next_probe = std::max<long long>(probe + 1LL, 2LL * probe + 1LL);
            if (next_probe > max_k) next_probe = max_k;
            probe = static_cast<int>(next_probe);
        }
    }

    if (lo > hi) {
        // Cached bounds can be tighter than the fresh lower bound
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
    if (lo - 1 > bounds.max_fail) {
        bounds.max_fail = lo - 1;
    }
    return lo;
}
