// Partner generation and diagonal index helpers
#include "partners.h"
#include "utils.h"
#include "profile.h"
#include "../flipdist/treedist.h" // for hasParentChildEdge
#include "../rotation_tree.h"
#include <algorithm>

using namespace flipdist;
#include <optional>
#include <array>

// Internal helper used only in this translation unit.
static std::vector<std::pair<int, int>> getIncidentEdgesLocal(const VectorRangeTreeMap &T, int node);
struct VertexEdgeCandidate
{
    std::pair<int,int> edge;
    int other_index;
};

static std::vector<VertexEdgeCandidate> collectVertexEdgeCandidates(const VectorRangeTreeMap &tree,
                                                                    int vertex)
{
    std::vector<VertexEdgeCandidate> result;
    auto edges = getIncidentEdgesLocal(tree, vertex);
    for (const auto &edge : edges)
    {
        int other = (edge.first == vertex) ? edge.second : edge.first;
        auto range = tree.getRange(other);
        int other_index = range.first;
        result.push_back({edge, other_index});
    }
    return result;
}

static bool isInternalDiagonal(const std::pair<int,int> &diag)
{
    return diag.second - diag.first > 1;
}

static DiagonalEdge makeDiagonalEdge(const VectorRangeTreeMap &T, int parent, int child)
{
    auto diag = T.getRange(child);
    if (diag.second <= diag.first)
    {
        int pos = diag.first;
        diag = {pos, pos + 1}; // boundary placeholder
    }
    return DiagonalEdge{diag, parent, child};
}

static std::vector<std::pair<int, int>> getIncidentEdgesLocal(const VectorRangeTreeMap &T, int node)
{
    std::vector<std::pair<int, int>> incident;

    try
    {
        if (!T.isOriginal(node))
            return incident;

        int left = T.getLeftChild(node);
        int right = T.getRightChild(node);
        if (left >= 0 && T.isOriginal(left))
        {
            incident.emplace_back(node, left);
        }
        if (right >= 0 && T.isOriginal(right))
        {
            incident.emplace_back(node, right);
        }

        int parent = T.getParent(node);
        if (parent >= 0 && T.isOriginal(parent))
        {
            incident.emplace_back(parent, node);
        }
    }
    catch (...)
    {
        // Return empty on error
    }

    return incident;
}

std::unordered_map<std::pair<int,int>, int, PairHash, PairEq>
buildEndpointIndex(const VectorRangeTreeMap &tree)
{
    std::unordered_map<std::pair<int,int>, int, PairHash, PairEq> index;
    index.reserve(tree.original_nodes.size());
    std::pair<int,int> fullRange{-1, -1};
    if (tree.root != VectorRangeTreeMap::NO_CHILD && tree.isOriginal(tree.root))
        fullRange = tree.getRange(tree.root);
    for (int node : tree.original_nodes)
    {
        auto endpoints = tree.diagonalEndpoints(node);
        if (endpoints.first < 0 || endpoints.second < 0) continue;
        if (endpoints.first >= endpoints.second) continue;
        // boundary edges (adjacent vertices) do not correspond to internal diagonals.
        if (endpoints.second - endpoints.first <= 1) continue;
        // The root's range corresponds to the outer boundary of the induced polygon,
        // not an internal diagonal.
        if (endpoints == fullRange) continue;
        index.emplace(endpoints, node);
    }
    return index;
}

std::unordered_map<std::pair<int,int>, int, PairHash, PairEq>
buildDiagonalNodeMap(const VectorRangeTreeMap &tree)
{
    return buildEndpointIndex(tree);
}

static std::string canonicalEdgePairKey(const std::pair<int,int> &a,
                                        const std::pair<int,int> &b)
{
    auto norm = [](std::pair<int,int> e) {
        if (e.first < e.second) return e;
        return std::make_pair(e.second, e.first);
    };
    auto ea = norm(a);
    auto eb = norm(b);
    if (ea > eb)
        std::swap(ea, eb);
    return std::to_string(ea.first) + "," + std::to_string(ea.second) + "|" +
           std::to_string(eb.first) + "," + std::to_string(eb.second);
}

static int polygonVertexCount(const VectorRangeTreeMap &tree)
{
    if (tree.original_nodes.empty())
        return 0;
    int maxVertex = 0;
    for (int node : tree.original_nodes)
    {
        auto range = tree.getRange(node);
        maxVertex = std::max({maxVertex, range.first, range.second});
    }
    return std::max(0, maxVertex + 1);
}

static std::optional<std::pair<std::pair<int,int>, std::pair<int,int>>>
pickBoundingEdges(const std::vector<VertexEdgeCandidate> &candidates,
                  int baseIndex,
                  int targetIndex,
                  int polygonVertices)
{
    if (candidates.size() < 2 || polygonVertices <= 0)
        return std::nullopt;

    auto normAngle = [&](int vertexIndex) -> int {
        long long raw = static_cast<long long>(vertexIndex) - static_cast<long long>(baseIndex);
        long long mod = raw % polygonVertices;
        if (mod < 0)
            mod += polygonVertices;
        return static_cast<int>(mod);
    };

    const int targetAngle = normAngle(targetIndex);

    std::vector<std::pair<int, const VertexEdgeCandidate *>> ordered;
    ordered.reserve(candidates.size());
    for (const auto &cand : candidates)
    {
        ordered.emplace_back(normAngle(cand.other_index), &cand);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    const VertexEdgeCandidate *lower = nullptr;
    const VertexEdgeCandidate *upper = nullptr;

    for (const auto &entry : ordered)
    {
        if (entry.first == targetAngle)
        {
            // Target diagonal already present; nothing to branch on.
            return std::nullopt;
        }
        if (entry.first < targetAngle)
        {
            lower = entry.second;
        }
        else if (entry.first > targetAngle && !upper)
        {
            upper = entry.second;
            break;
        }
    }

    if (!lower && !ordered.empty())
    {
        lower = ordered.back().second;
    }
    if (!upper && !ordered.empty())
    {
        upper = ordered.front().second;
    }

    if (!lower || !upper || lower->edge == upper->edge)
        return std::nullopt;

    return std::make_pair(lower->edge, upper->edge);
}

static std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>
enumerateAdjacentPairs(const std::vector<VertexEdgeCandidate> &candidates,
                       int baseIndex,
                       int polygonVertices)
{
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> result;
    if (candidates.size() < 2 || polygonVertices <= 0)
        return result;

    auto normAngle = [&](int vertexIndex) -> int {
        long long raw = static_cast<long long>(vertexIndex) - static_cast<long long>(baseIndex);
        long long mod = raw % polygonVertices;
        if (mod < 0)
            mod += polygonVertices;
        return static_cast<int>(mod);
    };

    std::vector<std::pair<int, const VertexEdgeCandidate *>> ordered;
    ordered.reserve(candidates.size());
    for (const auto &cand : candidates)
    {
        ordered.emplace_back(normAngle(cand.other_index), &cand);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    for (size_t i = 0; i < ordered.size(); ++i)
    {
        const auto *first = ordered[i].second;
        const auto *second = ordered[(i + 1) % ordered.size()].second;
        if (first->edge == second->edge)
            continue;
        result.emplace_back(first->edge, second->edge);
    }
    return result;
}

std::vector<std::pair<DiagonalEdge, DiagonalEdge>>
buildPartnerPairs(const VectorRangeTreeMap &start,
                  const VectorRangeTreeMap &target)
{
    flipdist::profile::ScopedTimer timer(&flipdist::profile::g.ns_buildPartnerPairs);
    if (flipdist::profile::enabled())
        ++flipdist::profile::g.calls_buildPartnerPairs;

    struct PartnerPair
    {
        DiagonalEdge first;
        DiagonalEdge second;
        int score;
        int insideA;
        int insideB;
    };
    std::vector<PartnerPair> partners;
    auto targetIndex = buildEndpointIndex(target);
    auto startIndex  = buildEndpointIndex(start);
    std::unordered_set<std::string> seen;
    const int polygonVertices = polygonVertexCount(start);

    auto orientEdge = [&](const std::pair<int,int> &edge)
        -> std::optional<DiagonalEdge>
    {
        if (hasParentChildEdge(start, edge.first, edge.second))
            return makeDiagonalEdge(start, edge.first, edge.second);
        if (hasParentChildEdge(start, edge.second, edge.first))
            return makeDiagonalEdge(start, edge.second, edge.first);
        return std::nullopt;
    };

    auto countInside = [&](const DiagonalEdge &edge) {
        int inside = 0;
        for (const auto &entry : targetIndex)
        {
            if (startIndex.count(entry.first))
                continue;
            if (entry.first.first >= edge.diag.first && entry.first.second <= edge.diag.second)
                ++inside;
        }
        return inside;
    };

    auto computeScore = [&](const DiagonalEdge &a, const DiagonalEdge &b,
                            int insideA, int insideB) {
        int widthA = a.diag.second - a.diag.first;
        int widthB = b.diag.second - b.diag.first;
        int span = widthA + widthB;
        int inside = insideA + insideB;
        return inside * 1000 + span;
    };

    auto addPair = [&](const DiagonalEdge &a, const DiagonalEdge &b) {
        if (!isInternalDiagonal(a.diag) || !isInternalDiagonal(b.diag))
            return;
        std::string key = canonicalEdgePairKey({a.parent,a.child}, {b.parent,b.child});
        if (!seen.insert(key).second)
            return;
        int insideA = countInside(a);
        int insideB = countInside(b);
        if (insideA == 0 && insideB == 0)
            return;
        int score = computeScore(a, b, insideA, insideB);
        partners.push_back({a, b, score, insideA, insideB});
    };

    auto emitCyclicPairs = [&](const std::vector<VertexEdgeCandidate> &infos, int baseVertex) {
        auto cyclic = enumerateAdjacentPairs(infos, baseVertex, polygonVertices);
        for (const auto &edges : cyclic)
        {
            auto maybeA = orientEdge(edges.first);
            auto maybeB = orientEdge(edges.second);
            if (maybeA && maybeB)
                addPair(*maybeA, *maybeB);
        }
    };

    // Limit to diagonals that are actual conflicts: present in start but missing in target.
    for (int node : start.original_nodes)
    {
        auto diag = start.diagonalEndpoints(node);
        if (diag.first < 0 || diag.second <= diag.first)
            continue;
        if (targetIndex.count(diag))
            continue; // only missing diagonals

        int left = diag.first;
        int right = diag.second;

        auto leftInfos = collectVertexEdgeCandidates(start, left);
        emitCyclicPairs(leftInfos, left);

        auto rightInfos = collectVertexEdgeCandidates(start, right);
        emitCyclicPairs(rightInfos, right);
    }

    for (const auto &entry : targetIndex)
    {
        const auto &diag = entry.first;
        if (startIndex.count(diag))
            continue; // only missing in start

        int left = diag.first;
        int right = diag.second;

        int leftIndex = left;
        int rightIndex = right;

        auto leftInfos = collectVertexEdgeCandidates(start, left);
        if (auto pair = pickBoundingEdges(leftInfos, leftIndex, rightIndex, polygonVertices))
        {
            auto maybeA = orientEdge(pair->first);
            auto maybeB = orientEdge(pair->second);
            if (maybeA && maybeB)
            {
                addPair(*maybeA, *maybeB);
            }
        }
        else if (leftInfos.size() >= 2)
        {
            for (size_t i = 0; i + 1 < leftInfos.size(); ++i)
            {
                for (size_t j = i + 1; j < leftInfos.size(); ++j)
                {
                    auto maybeA = orientEdge(leftInfos[i].edge);
                    auto maybeB = orientEdge(leftInfos[j].edge);
                    if (maybeA && maybeB)
                    {
                        addPair(*maybeA, *maybeB);
                    }
                }
            }
        }
        else if (leftInfos.size() == 1)
        {
            auto maybeA = orientEdge(leftInfos[0].edge);
            auto maybeB = orientEdge({left, right});
            if (maybeA && maybeB)
            {
                addPair(*maybeA, *maybeB);
            }
            else
            {
                int node = left;
                int parent = start.getParent(node);
                int sibling = VectorRangeTreeMap::NO_CHILD;
                if (parent != VectorRangeTreeMap::NO_PARENT && start.isOriginal(parent))
                {
                    int lch = start.getLeftChild(parent);
                    int rch = start.getRightChild(parent);
                    if (lch != node && lch != VectorRangeTreeMap::NO_CHILD) sibling = lch;
                    if (rch != node && rch != VectorRangeTreeMap::NO_CHILD) sibling = rch;
                    if (sibling != VectorRangeTreeMap::NO_CHILD && start.isOriginal(sibling))
                    {
                        auto maybeC = orientEdge({parent, sibling});
                        if (maybeA && maybeC)
                            addPair(*maybeA, *maybeC);
                    }
                }
            }
        }
        emitCyclicPairs(leftInfos, leftIndex);
        auto rightInfos = collectVertexEdgeCandidates(start, right);
        if (auto pair = pickBoundingEdges(rightInfos, rightIndex, leftIndex, polygonVertices))
        {
            auto maybeA = orientEdge(pair->first);
            auto maybeB = orientEdge(pair->second);
            if (maybeA && maybeB)
            {
                addPair(*maybeA, *maybeB);
            }
        }
        else if (rightInfos.size() >= 2)
        {
            for (size_t i = 0; i + 1 < rightInfos.size(); ++i)
            {
                for (size_t j = i + 1; j < rightInfos.size(); ++j)
                {
                    auto maybeA = orientEdge(rightInfos[i].edge);
                    auto maybeB = orientEdge(rightInfos[j].edge);
                    if (maybeA && maybeB)
                    {
                        addPair(*maybeA, *maybeB);
                    }
                }
            }
        }
        else if (rightInfos.size() == 1)
        {
            auto maybeA = orientEdge(rightInfos[0].edge);
            auto maybeB = orientEdge({right, left});
            if (maybeA && maybeB)
            {
                addPair(*maybeA, *maybeB);
            }
            else
            {
                int node = right;
                int parent = start.getParent(node);
                int sibling = VectorRangeTreeMap::NO_CHILD;
                if (parent != VectorRangeTreeMap::NO_PARENT && start.isOriginal(parent))
                {
                    int lch = start.getLeftChild(parent);
                    int rch = start.getRightChild(parent);
                    if (lch != node && lch != VectorRangeTreeMap::NO_CHILD) sibling = lch;
                    if (rch != node && rch != VectorRangeTreeMap::NO_CHILD) sibling = rch;
                    if (sibling != VectorRangeTreeMap::NO_CHILD && start.isOriginal(sibling))
                    {
                        auto maybeC = orientEdge({parent, sibling});
                        if (maybeA && maybeC)
                            addPair(*maybeA, *maybeC);
                    }
                }
            }
        }
        emitCyclicPairs(rightInfos, rightIndex);
    }

    auto diagWidth = [](const DiagonalEdge &edge) -> int {
        return edge.diag.second - edge.diag.first;
    };
    std::sort(partners.begin(), partners.end(),
              [&](const auto &lhs, const auto &rhs) {
                  if (lhs.score != rhs.score)
                      return lhs.score > rhs.score;
                  int lw = std::max(diagWidth(lhs.first), diagWidth(lhs.second));
                  int rw = std::max(diagWidth(rhs.first), diagWidth(rhs.second));
                  if (lw != rw)
                      return lw > rw;
                  int lw2 = std::min(diagWidth(lhs.first), diagWidth(lhs.second));
                  int rw2 = std::min(diagWidth(rhs.first), diagWidth(rhs.second));
                  if (lw2 != rw2)
                      return lw2 > rw2;
                  if (lhs.first.diag != rhs.first.diag)
                      return lhs.first.diag < rhs.first.diag;
                  return lhs.second.diag < rhs.second.diag;
              });

    size_t missingDiagonals = 0;
    for (const auto &entry : targetIndex)
    {
        if (!startIndex.count(entry.first))
            ++missingDiagonals;
    }
    // Partner pairs act as a search accelerator, but overly aggressive caps can
    // make the Li–Xia branching incomplete (missing a required wedge).
    // Keep this bounded for performance, but allow enough headroom to cover
    // most missing diagonals in medium-size instances.
    // NOTE: Overly small caps have proven to break completeness on some n=12
    // parity cases against the Java BFS oracle (we can miss the needed wedge).
    // Keep a cap for performance, but give more headroom at moderate sizes.
    const size_t partnerCap = std::min<size_t>(128, std::max<size_t>(12, missingDiagonals * 6));

    std::unordered_map<std::pair<int,int>, size_t, PairHash, PairEq> perDiagCount;
    perDiagCount.reserve(partners.size());
    const size_t perDiagLimit = 4;
    std::unordered_set<std::string> diagSeen;
    diagSeen.reserve(partners.size());
    std::vector<std::pair<DiagonalEdge, DiagonalEdge>> ordered;
    ordered.reserve(std::min(partners.size(), partnerCap));
    std::vector<PartnerPair> fallbackPairs;
    fallbackPairs.reserve(partners.size());
    for (const auto &pair : partners)
    {
        std::string key = makeDiagonalPairKey(pair.first, pair.second);
        if (!diagSeen.insert(key).second)
            continue;
        auto bump = [&](const DiagonalEdge &edge) -> bool {
            auto &count = perDiagCount[edge.diag];
            if (count >= perDiagLimit)
                return false;
            ++count;
            return true;
        };
        if (!bump(pair.first) || !bump(pair.second))
            continue;
        if (pair.insideA > 0 && pair.insideB > 0)
        {
            ordered.emplace_back(pair.first, pair.second);
            if (ordered.size() >= partnerCap)
                break;
        }
        else
        {
            fallbackPairs.push_back(pair);
        }
    }
    if (ordered.size() < partnerCap && !fallbackPairs.empty())
    {
        size_t keep = std::min(partnerCap - ordered.size(), fallbackPairs.size());
        std::partial_sort(fallbackPairs.begin(), fallbackPairs.begin() + keep, fallbackPairs.end(),
                          [](const PartnerPair &a, const PartnerPair &b) {
                              int ga = a.insideA + a.insideB;
                              int gb = b.insideA + b.insideB;
                              if (ga != gb) return ga > gb;
                              return a.score > b.score;
                          });
        for (size_t i = 0; i < keep; ++i)
        {
            ordered.emplace_back(fallbackPairs[i].first, fallbackPairs[i].second);
        }
    }

    if (flipdist::profile::enabled())
        flipdist::profile::g.total_partnerPairs += static_cast<std::uint64_t>(ordered.size());
    return ordered;
}
