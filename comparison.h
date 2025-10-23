// Declares the data structures and helpers used to describe a single solver
// run. `ComparisonOptions` captures CLI-style switches, `ComparisonRow` is the
// serialisable telemetry record, and `runComparison`/`runBidirectional`
// implement the distance computation with optional bidirectional fallback.

#ifndef COMPARISON_H
#define COMPARISON_H

#include "rotation_tree.h"
#include <cstddef>
#include <cstdint>
#include <string>

struct ComparisonOptions {
    std::string program = "test_asan";
    std::string case_type = "unknown";
    std::string direction = "a->b";
    long long seed = -1;
    double time_limit_sec = 5.0;
    std::size_t visited_cap = 5'000'000;
    std::size_t queue_cap = 5'000'000;
    bool use_bidir_on_timeout = false;
    std::size_t bidir_state_cap = 5'000'000;
    bool prefer_bidir = false;
};

struct ComparisonRow {
    std::string program;
    std::string case_type;
    int n;
    long long seed;
    std::string direction;
    int distance;
    double time_ms;
    std::size_t expanded;
    std::size_t enqueued;
    std::size_t visited;
    std::size_t max_queue;
    std::size_t duplicates;
    std::string status;
    std::string tree_a;
    std::string tree_b;
    std::string solver;
};

ComparisonRow runComparison(const VectorRangeTreeMap& start,
                            const VectorRangeTreeMap& target,
                            const ComparisonOptions& opts);

std::string comparisonRowToJson(const ComparisonRow& row);

struct ComparisonPair {
    ComparisonRow forward;
    ComparisonRow reverse;
    bool distance_agrees;
};

ComparisonPair runBidirectional(const VectorRangeTreeMap& a,
                                const VectorRangeTreeMap& b,
                                const ComparisonOptions& opts);

#endif // COMPARISON_H
