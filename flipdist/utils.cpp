// Utility functions used across FlipDist
#include "utils.h"
#include "../rotation_tree.h" // for PairHash/PairEq definitions
#include <algorithm>

std::string canonicalEdgeListKey(const std::vector<std::pair<int,int>> &edges)
{
    if (edges.empty())
        return {};
    std::vector<std::pair<int,int>> sorted = edges;
    std::sort(sorted.begin(), sorted.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.first != rhs.first) return lhs.first < rhs.first;
        return lhs.second < rhs.second;
    });
    std::string key;
    key.reserve(sorted.size() * 8);
    bool first = true;
    for (const auto &edge : sorted)
    {
        if (!first)
            key.push_back(';');
        first = false;
        key += std::to_string(edge.first);
        key.push_back(',');
        key += std::to_string(edge.second);
    }
    return key;
}

std::string canonicalDiagonalPairListKey(const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &pairs)
{
    if (pairs.empty())
        return {};

    auto encodeEdge = [](const DiagonalEdge &edge) {
        return std::array<int, 4>{edge.diag.first, edge.diag.second, edge.parent, edge.child};
    };

    std::vector<std::array<int, 8>> normalized;
    normalized.reserve(pairs.size());
    for (const auto &pair : pairs)
    {
        auto first  = encodeEdge(pair.first);
        auto second = encodeEdge(pair.second);
        if (second < first)
            std::swap(first, second);
        normalized.push_back({first[0], first[1], first[2], first[3],
                              second[0], second[1], second[2], second[3]});
    }

    std::sort(normalized.begin(), normalized.end());

    std::string key;
    key.reserve(normalized.size() * 32);
    bool firstEntry = true;
    for (const auto &entry : normalized)
    {
        if (!firstEntry)
            key.push_back(';');
        firstEntry = false;
        for (size_t i = 0; i < entry.size(); ++i)
        {
            if (i > 0)
                key.push_back(',');
            key += std::to_string(entry[i]);
        }
    }
    return key;
}

std::string canonicalDiagonalSetKey(const std::unordered_set<std::pair<int,int>, PairHash, PairEq> &diags)
{
    if (diags.empty())
        return {};
    std::vector<std::pair<int,int>> sorted(diags.begin(), diags.end());
    std::sort(sorted.begin(), sorted.end());
    std::string key;
    key.reserve(sorted.size() * 8);
    bool first = true;
    for (const auto &d : sorted)
    {
        if (!first)
            key.push_back(';');
        first = false;
        key += std::to_string(d.first);
        key.push_back(',');
        key += std::to_string(d.second);
    }
    return key;
}

std::string makeDiagonalPairKey(const DiagonalEdge &a,
                                const DiagonalEdge &b)
{
    auto diagA = a.diag;
    auto diagB = b.diag;
    if (diagB < diagA)
        std::swap(diagA, diagB);
    return std::to_string(diagA.first) + "," + std::to_string(diagA.second) + "|" +
           std::to_string(diagB.first) + "," + std::to_string(diagB.second);
}

std::string makeRangeMemoKey(const std::string &treePairSig,
                             const std::string &conflictSig,
                             const std::pair<int,int> &range)
{
    std::string key = treePairSig;
    key += "|R:";
    key += std::to_string(range.first);
    key.push_back(',');
    key += std::to_string(range.second);
    key += "|C:";
    key += conflictSig;
    return key;
}
