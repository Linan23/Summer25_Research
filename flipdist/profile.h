#pragma once

#include <cstdint>
#include <chrono>
#include <ostream>

namespace flipdist::profile {

struct Stats
{
    std::uint64_t calls_flipDistTree{0};
    std::uint64_t calls_treeDistI{0};
    std::uint64_t calls_treeDistS{0};
    std::uint64_t calls_treeSignature{0};
    std::uint64_t calls_collectConflicts{0};
    std::uint64_t calls_buildPartnerPairs{0};
    std::uint64_t calls_findFreeEdges{0};
    std::uint64_t calls_safeCopyTree{0};
    std::uint64_t calls_partitionAlongEdge{0};
    std::uint64_t calls_minRotationsBFS_internal{0};
    std::uint64_t memo_hits_minRotationsBFS_internal{0};
    std::uint64_t successes_minRotationsBFS_internal{0};
    std::uint64_t failures_minRotationsBFS_internal{0};

    std::uint64_t memo_hits_flipDistTree{0};
    std::uint64_t memo_hits_treeDistI{0};
    std::uint64_t memo_hits_treeDistS{0};

    std::uint64_t total_partnerPairs{0};
    std::uint64_t total_freeEdgeCandidates{0};
    std::uint64_t total_freeEdgeAttempts{0};
    std::uint64_t total_partitionMismatches{0};

    std::uint64_t max_conflicts_seen{0};

    std::uint64_t ns_flipDistTree{0};
    std::uint64_t ns_treeDistI{0};
    std::uint64_t ns_treeDistS{0};
    std::uint64_t ns_treeSignature{0};
    std::uint64_t ns_collectConflicts{0};
    std::uint64_t ns_buildPartnerPairs{0};
    std::uint64_t ns_findFreeEdges{0};
    std::uint64_t ns_safeCopyTree{0};
    std::uint64_t ns_partitionAlongEdge{0};
    std::uint64_t ns_minRotationsBFS_internal{0};
};

extern Stats g;

bool enabled();
void reset();
void print(std::ostream &os, const char *tag = nullptr);
void heartbeat(const char *tag = nullptr);

struct ScopedTimer
{
    std::uint64_t *dest_ns;
    std::chrono::steady_clock::time_point start;
    bool active;

    explicit ScopedTimer(std::uint64_t *dest);
    ~ScopedTimer();
};

} // namespace flipdist::profile
