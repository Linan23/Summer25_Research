#pragma once

#include "bf_bst.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

extern const bool DEBUG;

struct ProfileStats {
    bool enabled = false;
    long long calls_flipdist = 0;
    long long calls_tdi = 0;
    long long calls_tds = 0;
    long long calls_partition = 0;
    long long free_edge_hits = 0;
    long long free_edge_misses = 0;
    long long s_branch_calls = 0;
    long long s_empty_branch_calls = 0;
    long long independent_subsets_initial = 0;
    long long independent_subsets_s = 0;
    long long max_s_size = 0;
    long long max_i_size = 0;
    double time_flipdist_ms = 0.0;
    double time_tdi_ms = 0.0;
    double time_tds_ms = 0.0;
    double time_partition_ms = 0.0;
    std::chrono::steady_clock::time_point start_time;
    int abort_ms = -1;
};

extern ProfileStats g_profile;

struct Key128 {
    std::uint64_t hi = 0;
    std::uint64_t lo = 0;
    bool operator==(const Key128 &o) const {
        return hi == o.hi && lo == o.lo;
    }
};

struct Key128Hash {
    std::size_t operator()(const Key128 &k) const;
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

void dedupeSPairs(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &pairs);

Key128 makeKeyPair(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);
Key128 makeKeyFlip(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k);
Key128 makeKeyI(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final, int k,
                const std::vector<std::pair<int, int>> &I);
Key128 makeKeyS(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end, int k,
                const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);
Key128 makeKeySBase(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_end,
                    const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &S);
Key128 makeKeyIBase(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final,
                    const std::vector<std::pair<int, int>> &I);

bool tryBoundsPrune(const std::unordered_map<Key128, KBounds, Key128Hash> &bounds_map,
                    const Key128 &key, int k, bool &value_out);
void updateBounds(std::unordered_map<Key128, KBounds, Key128Hash> &bounds_map,
                  const Key128 &key, int k, bool value);
int requiredBudgetFromBounds(const std::unordered_map<Key128, KBounds, Key128Hash> &bounds_map,
                             const Key128 &key);

void resetMemo();
void initProfile();
bool profileAbortRequested();

struct ScopedTimer {
    double *acc = nullptr;
    std::chrono::steady_clock::time_point t0;
    explicit ScopedTimer(double *a);
    ~ScopedTimer();
};

void debugPrint(const std::string &msg);
