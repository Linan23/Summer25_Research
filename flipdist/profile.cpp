#include "profile.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace flipdist::profile {

Stats g;

bool enabled()
{
    static int cached = -1;
    if (cached == -1)
        cached = (std::getenv("FLIPDIST_PROFILE") != nullptr) ? 1 : 0;
    return cached == 1;
}

void reset()
{
    g = Stats{};
}

static double nsToMs(std::uint64_t ns)
{
    return static_cast<double>(ns) / 1'000'000.0;
}

void print(std::ostream &os, const char *tag)
{
    if (!enabled())
        return;

    if (tag && *tag)
        os << "\n[PROFILE] " << tag << "\n";
    else
        os << "\n[PROFILE]\n";

    auto line = [&](const char *name, std::uint64_t calls, std::uint64_t ns) {
        double ms = nsToMs(ns);
        double avg_us = (calls > 0) ? (ms * 1000.0 / static_cast<double>(calls)) : 0.0;
        os << "  - " << std::left << std::setw(24) << name
           << " calls=" << std::right << std::setw(10) << calls
           << " time_ms=" << std::setw(12) << std::fixed << std::setprecision(3) << ms
           << " avg_us=" << std::setw(10) << std::fixed << std::setprecision(2) << avg_us
           << "\n";
    };

    line("FlipDistTree", g.calls_flipDistTree, g.ns_flipDistTree);
    line("TreeDistI", g.calls_treeDistI, g.ns_treeDistI);
    line("TreeDistS", g.calls_treeDistS, g.ns_treeDistS);
    line("treeSignature", g.calls_treeSignature, g.ns_treeSignature);
    line("collectConflicts", g.calls_collectConflicts, g.ns_collectConflicts);
    line("buildPartnerPairs", g.calls_buildPartnerPairs, g.ns_buildPartnerPairs);
    line("findFreeEdges", g.calls_findFreeEdges, g.ns_findFreeEdges);
    line("safeCopyTree", g.calls_safeCopyTree, g.ns_safeCopyTree);
    line("partitionAlongEdge", g.calls_partitionAlongEdge, g.ns_partitionAlongEdge);
    line("MinRotationsBFS", g.calls_minRotationsBFS_internal, g.ns_minRotationsBFS_internal);
    if (g.calls_minRotationsBFS_internal > 0)
    {
        os << "    * bfs_memo_hits=" << g.memo_hits_minRotationsBFS_internal
           << " bfs_ok=" << g.successes_minRotationsBFS_internal
           << " bfs_fail=" << g.failures_minRotationsBFS_internal
           << "\n";
    }

    os << "  - memo_hits FlipDistTree=" << g.memo_hits_flipDistTree
       << " TreeDistI=" << g.memo_hits_treeDistI
       << " TreeDistS=" << g.memo_hits_treeDistS
       << "\n";
    os << "  - totals partnerPairs=" << g.total_partnerPairs
       << " freeEdgeCandidates=" << g.total_freeEdgeCandidates
       << " freeEdgeAttempts=" << g.total_freeEdgeAttempts
       << " partitionMismatches=" << g.total_partitionMismatches
       << "\n";
    os << "  - max_conflicts_seen=" << g.max_conflicts_seen << "\n";
}

void heartbeat(const char *tag)
{
    if (!enabled())
        return;

    static auto last = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    if (now - last < std::chrono::seconds(5))
        return;
    last = now;
    print(std::cerr, tag);
}

ScopedTimer::ScopedTimer(std::uint64_t *dest)
    : dest_ns(dest), start(std::chrono::steady_clock::now()), active(enabled() && dest != nullptr)
{
}

ScopedTimer::~ScopedTimer()
{
    if (!active)
        return;
    auto end = std::chrono::steady_clock::now();
    *dest_ns += static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
}

} // namespace flipdist::profile
