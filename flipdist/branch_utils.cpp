// Branching helpers for FlipDist
#include "branch_utils.h"
#include "utils.h"
#include <unordered_set>
#include <set>

namespace {
constexpr bool PROFILE = false; // profiling disabled by default

struct DedupProfileEntry
{
    long long calls = 0;
    long long before = 0;
    long long after = 0;
};

struct BranchProfileEntry
{
    long long calls = 0;
    long long size_before = 0;
    long long size_after = 0;
};

static std::unordered_map<std::string, DedupProfileEntry> g_dedupProfiles;
static std::unordered_map<std::string, BranchProfileEntry> g_branchProfiles;

static inline bool isInternalDiagonal(const std::pair<int,int> &diag)
{
    return diag.second - diag.first > 1;
}

static void recordDedupSample(const char *tag, size_t before, size_t after)
{
    if (!PROFILE || tag == nullptr)
        return;
    auto &entry = g_dedupProfiles[tag];
    entry.calls += 1;
    entry.before += static_cast<long long>(before);
    entry.after += static_cast<long long>(after);
}
} // namespace

bool appendValidPair(std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &dest,
                     const DiagonalEdge &a,
                     const DiagonalEdge &b)
{
    const bool aInternal = isInternalDiagonal(a.diag);
    const bool bInternal = isInternalDiagonal(b.diag);
    // Li–Xia's S-pairs can include polygon boundary edges; in this binary-tree
    // model those show up as width-1 diagonals / placeholder edges that are
    // not rotatable. Keep pairs with at least one internal diagonal so the
    // branching still considers the rotatable edge.
    if (!aInternal && !bInternal)
        return false;
    if (a.diag == b.diag &&
        ((a.parent == b.parent && a.child == b.child) ||
         (a.parent == b.child && a.child == b.parent)))
        return false;
    dest.emplace_back(a, b);
    return true;
}

void deduplicateDiagonalPairs(std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &pairs,
                              const char *tag)
{
    size_t before = pairs.size();
    if (before <= 1)
    {
        recordDedupSample(tag, before, before);
        return;
    }

    std::unordered_set<std::string> seen;
    seen.reserve(pairs.size() * 2);
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> filtered;
    filtered.reserve(pairs.size());
    for (const auto &p : pairs)
    {
        std::string key = makeDiagonalPairKey(p.first, p.second);
        if (seen.insert(key).second)
            filtered.push_back(p);
    }
    if (filtered.size() != pairs.size())
        pairs.swap(filtered);
    recordDedupSample(tag, before, pairs.size());
}

void recordBranchSample(const char *tag, size_t before, size_t after)
{
    if (!PROFILE || tag == nullptr)
        return;
    auto &entry = g_branchProfiles[tag];
    entry.calls += 1;
    entry.size_before += static_cast<long long>(before);
    entry.size_after += static_cast<long long>(after);
}

std::pair<std::vector<std::pair<DiagonalEdge, DiagonalEdge>>,
          std::vector<std::pair<DiagonalEdge, DiagonalEdge>>>
partitionS(const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &S,
           const VectorRangeTreeMap &T1, const VectorRangeTreeMap &T2)
{
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> S1, S2;

    // Get node sets for each partition
    std::set<int> nodes1, nodes2;
    for (int node : T1.original_nodes)
        nodes1.insert(node);
    for (int node : T2.original_nodes)
        nodes2.insert(node);

    for (const auto &edgePair : S)
    {
        const auto &edge1 = edgePair.first;
        const auto &edge2 = edgePair.second;

        auto edge_in_nodes = [](const DiagonalEdge &edge, const std::set<int> &nodes) -> bool
        {
            // Boundary/placeholder edges are represented with one missing endpoint
            // (child < 0). Assign them based on the incident internal node so
            // S-pairs containing a boundary edge remain in the correct subtree.
            if (edge.child < 0)
                return nodes.count(edge.parent) != 0;
            if (edge.parent < 0)
                return nodes.count(edge.child) != 0;
            return nodes.count(edge.parent) && nodes.count(edge.child);
        };

        bool edge1_in_T1 = edge_in_nodes(edge1, nodes1);
        bool edge2_in_T1 = edge_in_nodes(edge2, nodes1);

        bool edge1_in_T2 = edge_in_nodes(edge1, nodes2);
        bool edge2_in_T2 = edge_in_nodes(edge2, nodes2);

        if (edge1_in_T1 && edge2_in_T1)
        {
            S1.push_back(edgePair);
        }
        else if (edge1_in_T2 && edge2_in_T2)
        {
            S2.push_back(edgePair);
        }
        // cross-partition pairs are ignored
    }

    deduplicateDiagonalPairs(S1, "partitionS_S1");
    deduplicateDiagonalPairs(S2, "partitionS_S2");

    return {S1, S2};
}
