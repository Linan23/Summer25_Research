// Precomputed all-pairs rotation distances for very small BSTs (n <= 9).
//
// The Li–Xia recursion frequently generates tiny subproblems. A naive BFS per
// query is cheap in isolation but becomes a major bottleneck when invoked tens
// of thousands of times. For n <= 9, the rotation graph has at most Catalan(9)
// = 4862 nodes, so we can precompute exact distances once and answer queries in
// O(1).

#include "small_bfs.h"

#include "conflicts.h"
#include "memo.h"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace flipdist {
namespace {

constexpr int kMaxSmallNodes = 9;

static std::string preorderKey(const std::vector<int> &pre)
{
    std::string key;
    key.reserve(pre.size() * 3);
    for (size_t i = 0; i < pre.size(); ++i)
    {
        if (i)
            key.push_back(',');
        key.append(std::to_string(pre[i]));
    }
    return key;
}

static std::vector<int> canonicalizePreorder(const std::vector<int> &pre,
                                             const std::vector<int> &in)
{
    std::vector<int> sorted = in;
    std::sort(sorted.begin(), sorted.end());
    std::vector<int> out;
    out.reserve(pre.size());
    for (int v : pre)
    {
        auto it = std::lower_bound(sorted.begin(), sorted.end(), v);
        if (it == sorted.end() || *it != v)
            return {};
        out.push_back(static_cast<int>(it - sorted.begin()) + 1);
    }
    return out;
}

static const std::vector<std::vector<int>> &generateAllBSTPreorders(int n)
{
    static std::vector<std::vector<std::vector<int>>> memo(kMaxSmallNodes + 1);
    static std::vector<char> built(kMaxSmallNodes + 1, 0);
    if (n < 0 || n > kMaxSmallNodes)
        return memo[0];
    if (built[n])
        return memo[n];

    std::vector<std::vector<int>> result;
    if (n == 0)
    {
        result.push_back({});
        memo[n] = result;
        built[n] = 1;
        return memo[n];
    }

    for (int root = 1; root <= n; ++root)
    {
        const auto &leftList = generateAllBSTPreorders(root - 1);
        const auto &rightList = generateAllBSTPreorders(n - root);

        for (const auto &left : leftList)
        {
            for (const auto &right : rightList)
            {
                std::vector<int> pre;
                pre.reserve(n);
                pre.push_back(root);
                pre.insert(pre.end(), left.begin(), left.end());
                for (int v : right)
                    pre.push_back(v + root);
                result.push_back(std::move(pre));
            }
        }
    }

    memo[n] = std::move(result);
    built[n] = 1;
    return memo[n];
}

struct SmallDistanceTable
{
    int n = 0;
    int count = 0;
    std::unordered_map<std::string, int> indexByPre;
    std::vector<std::vector<int>> adj;
    std::vector<std::uint8_t> dist; // flattened count*count matrix; 255=unreachable
};

static SmallDistanceTable buildTable(int n)
{
    SmallDistanceTable table;
    table.n = n;

    const auto &allPre = generateAllBSTPreorders(n);
    table.count = static_cast<int>(allPre.size());
    table.indexByPre.reserve(allPre.size() * 2);

    for (int i = 0; i < table.count; ++i)
    {
        table.indexByPre.emplace(preorderKey(allPre[i]), i);
    }

    table.adj.assign(table.count, {});

    std::vector<int> inorder(n);
    for (int i = 0; i < n; ++i)
        inorder[i] = i + 1;

    for (int i = 0; i < table.count; ++i)
    {
        VectorRangeTreeMap T;
        T.build(allPre[i], inorder);
        auto edges = getInternalEdges(T);
        auto &nbrs = table.adj[i];
        nbrs.reserve(edges.size());
        for (auto [p, c] : edges)
        {
            VectorRangeTreeMap nx = T;
            if (nx.getLeftChild(p) == c)
                nx.rotateRight(p);
            else if (nx.getRightChild(p) == c)
                nx.rotateLeft(p);
            else
                continue;

            auto [preNx, inNx] = serializeTree(nx);
            if (static_cast<int>(preNx.size()) != n || static_cast<int>(inNx.size()) != n)
                continue;
            auto it = table.indexByPre.find(preorderKey(preNx));
            if (it == table.indexByPre.end())
                continue;
            nbrs.push_back(it->second);
        }
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    const int N = table.count;
    table.dist.assign(static_cast<size_t>(N) * static_cast<size_t>(N), static_cast<std::uint8_t>(255));
    std::deque<int> q;
    std::vector<std::uint8_t> d(static_cast<size_t>(N), static_cast<std::uint8_t>(255));

    for (int src = 0; src < N; ++src)
    {
        std::fill(d.begin(), d.end(), static_cast<std::uint8_t>(255));
        q.clear();
        d[src] = 0;
        q.push_back(src);
        while (!q.empty())
        {
            int u = q.front();
            q.pop_front();
            std::uint8_t du = d[u];
            for (int v : table.adj[u])
            {
                if (d[v] != 255)
                    continue;
                d[v] = static_cast<std::uint8_t>(du + 1);
                q.push_back(v);
            }
        }
        for (int dst = 0; dst < N; ++dst)
        {
            table.dist[static_cast<size_t>(src) * static_cast<size_t>(N) + static_cast<size_t>(dst)] = d[dst];
        }
    }

    return table;
}

static const SmallDistanceTable *getTableForSize(int n)
{
    if (n < 0 || n > kMaxSmallNodes)
        return nullptr;
    static std::unordered_map<int, SmallDistanceTable> tables;
    auto it = tables.find(n);
    if (it != tables.end())
        return &it->second;
    tables.emplace(n, buildTable(n));
    return &tables.find(n)->second;
}

} // namespace

int smallRotationDistance(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap)
{
    if (cap < 0)
        return -2;
    auto [preA, inA] = serializeTree(A);
    auto [preB, inB] = serializeTree(B);
    if (preA.size() != preB.size())
        return -2;
    const int n = static_cast<int>(preA.size());
    if (n == 0)
        return 0;
    if (n > kMaxSmallNodes)
        return -2;

    auto canonA = canonicalizePreorder(preA, inA);
    auto canonB = canonicalizePreorder(preB, inB);
    if (canonA.empty() || canonB.empty())
        return -2;

    const auto *table = getTableForSize(n);
    if (!table)
        return -2;

    auto itA = table->indexByPre.find(preorderKey(canonA));
    auto itB = table->indexByPre.find(preorderKey(canonB));
    if (itA == table->indexByPre.end() || itB == table->indexByPre.end())
        return -2;

    const int idxA = itA->second;
    const int idxB = itB->second;
    const int N = table->count;
    const std::uint8_t d = table->dist[static_cast<size_t>(idxA) * static_cast<size_t>(N) + static_cast<size_t>(idxB)];
    if (d == 255)
        return -2;
    return (d <= cap) ? static_cast<int>(d) : -1;
}

} // namespace flipdist
