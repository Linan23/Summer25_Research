#include "memoization.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>

// Debug flag
const bool DEBUG = false;

ProfileStats g_profile;
std::unordered_map<Key128, bool, Key128Hash> g_memo_flipdist;
std::unordered_map<Key128, bool, Key128Hash> g_memo_tdi;
std::unordered_map<Key128, bool, Key128Hash> g_memo_tds;
std::unordered_map<Key128, KBounds, Key128Hash> g_kbounds;
std::unordered_map<Key128, KBounds, Key128Hash> g_tds_bounds;
std::unordered_map<Key128, KBounds, Key128Hash> g_tdi_bounds;

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
        uint64_t k1 = edgeKey(p.first);
        uint64_t k2 = edgeKey(p.second);
        if (k2 < k1) std::swap(k1, k2);
        Key128 key{k1, k2};
        if (seen.insert(key).second) {
            out.push_back(p);
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
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());
    TreeFingerprint out{1469598103934665603ULL, 1099511628211ULL};
    mixFingerprint(out, static_cast<uint64_t>(nodes.size()));
    mixFingerprint(out, static_cast<uint64_t>(T.root + 1));
    for (int v : nodes) {
        auto r = T.getRange(v);
        int L = T.getLeftChild(v);
        int R = T.getRightChild(v);
        int P = T.getParent(v);
        mixFingerprint(out, static_cast<uint64_t>(v + 1));
        mixFingerprint(out, (static_cast<uint64_t>(static_cast<uint32_t>(r.first)) << 32) ^
                            static_cast<uint32_t>(r.second));
        mixFingerprint(out, (static_cast<uint64_t>(static_cast<uint32_t>(L + 2)) << 32) ^
                            static_cast<uint32_t>(R + 2));
        mixFingerprint(out, static_cast<uint64_t>(P + 2));
    }
    return out;
}

Key128 makeKeyPair(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final) {
    TreeFingerprint a = treeFingerprint(T_init);
    TreeFingerprint b = treeFingerprint(T_final);
    Key128 key{1469598103934665603ULL, 1099511628211ULL};
    mixKey(key, 1);
    mixKey(key, a.h1); mixKey(key, a.h2);
    mixKey(key, b.h1); mixKey(key, b.h2);
    return key;
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
    if (value) {
        if (b.min_success < 0 || k < b.min_success) b.min_success = k;
        if (b.max_fail >= b.min_success) b.max_fail = b.min_success - 1;
    } else {
        if (k > b.max_fail) b.max_fail = k;
    }
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
}

void initProfile() {
    if (g_profile.enabled) return;
    const char *env = std::getenv("FLIPDIST_PROFILE");
    if (env && std::string(env) == "1") {
        g_profile.enabled = true;
        g_profile.start_time = std::chrono::steady_clock::now();
        const char *abort_env = std::getenv("FLIPDIST_PROFILE_ABORT_MS");
        if (abort_env) {
            g_profile.abort_ms = std::atoi(abort_env);
        }
    }
}

bool profileAbortRequested() {
    if (!g_profile.enabled || g_profile.abort_ms <= 0) return false;
    auto now = std::chrono::steady_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(now - g_profile.start_time).count();
    return elapsed_ms > g_profile.abort_ms;
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
