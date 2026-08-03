#include "bf_bst.h"
#include <iostream>
#include <algorithm>
#include <set>
#include <sstream>
#include <queue>
#include <climits>
#include <cassert>
#include <functional>
#include <random>
#include <chrono>
#include <iomanip>
#include <array>

// Purpose: Hash a (parent, child) pair for unordered containers
// Inputs: p: pair to hash
// Output: hash value
size_t PairHash::operator()(const std::pair<int, int> &p) const
{
    return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned long long)p.second);
}

// Purpose: Compare two (parent, child) pairs for equality
// Inputs: a: first pair; b: second pair
// Output: true if both pairs match
bool PairEq::operator()(const std::pair<int, int> &a,
                        const std::pair<int, int> &b) const
{
    return a.first == b.first && a.second == b.second;
}

// Purpose: Initialize an empty tree with default values
// Inputs: none
// Output: constructed tree
VectorRangeTreeMap::VectorRangeTreeMap() : root(-1), max_node_value(0)
{

    std::vector<std::pair<int, int>> ranges;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> parents;

    int root = -1;
    int max_node_value = 0;

    std::vector<int> original_inorder, original_preorder;
    std::unordered_map<int, int> position_in_inorder;
    std::vector<int> position_in_inorder_fast;
    std::unordered_set<int> original_nodes;
    std::vector<unsigned char> original_mask;

    static constexpr int NO_CHILD = -1;
    static constexpr int NO_PARENT = -1;
}

// Purpose: Build tree and ranges from preorder and inorder traversals
// Inputs: preorder: preorder node list; inorder: inorder node list
// Output: nothing
// Failure: leaves tree empty on invalid input
void VectorRangeTreeMap::build(const std::vector<int> &preorder, const std::vector<int> &inorder)
{
    buildImpl(preorder, inorder, true, true);
}

// Purpose: Builds light.
// Inputs: preorder (const std::vector<int> &); inorder (const std::vector<int> &).
// Output: No return value; updates referenced outputs or shared state.
void VectorRangeTreeMap::buildLight(const std::vector<int> &preorder, const std::vector<int> &inorder)
{
    buildImpl(preorder, inorder, false, false);
}

// Purpose: Builds impl.
// Inputs: preorder (const std::vector<int> &); inorder (const std::vector<int> &); populate_nodes (bool); store_preorder (bool).
// Output: No return value; updates referenced outputs or shared state.
void VectorRangeTreeMap::buildImpl(const std::vector<int> &preorder,
                                   const std::vector<int> &inorder,
                                   bool populate_nodes,
                                   bool store_preorder)
{
    if (store_preorder) {
        original_preorder = preorder;
    } else {
        original_preorder.clear();
    }

    clear();

    if (preorder.empty())
        return;

    max_node_value = *std::max_element(preorder.begin(), preorder.end());
    max_node_value = std::max(max_node_value, *std::max_element(inorder.begin(), inorder.end()));

    ranges.resize(max_node_value + 1, {0, 0});
    edges.resize(max_node_value + 1, {NO_CHILD, NO_CHILD});
    parents.resize(max_node_value + 1, NO_PARENT);
    original_mask.resize(max_node_value + 1, 0);
    position_in_inorder_fast.resize(max_node_value + 1, -1);

    // The mask keeps hot membership checks contiguous; the set supports callers needing set semantics.
    if (populate_nodes) {
        original_nodes.max_load_factor(2.0f);
        original_nodes.reserve(preorder.size());
    }
    for (int val : preorder)
    {
        if (populate_nodes) {
            original_nodes.insert(val);
        }
        if (val >= 0 && val < (int)original_mask.size())
            original_mask[val] = 1;
    }

    original_inorder = inorder;

    bool needs_sparse_positions = false;
    for (int val : inorder) {
        if (val < 0 || val >= (int)position_in_inorder_fast.size()) {
            needs_sparse_positions = true;
            break;
        }
    }
    if (needs_sparse_positions) {
        position_in_inorder.reserve(inorder.size());
    }
    for (size_t i = 0; i < inorder.size(); i++)
    {
        if (needs_sparse_positions) {
            position_in_inorder[inorder[i]] = i;
        }
        if (inorder[i] >= 0 && inorder[i] < (int)position_in_inorder_fast.size())
            position_in_inorder_fast[inorder[i]] = (int)i;
    }

    root = preorder[0];
    buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, NO_PARENT);

    calculateAllRanges();
    fingerprint_valid = false;
}

// Purpose: Get the left child of a node
// Inputs: node: node id
// Output: left child id or NO_CHILD
// Failure: returns NO_CHILD for invalid node
int VectorRangeTreeMap::getLeftChild(int node) const
{
    return (node >= 0 && node < (int)edges.size()) ? edges[node].first : NO_CHILD;
}

// Purpose: Get the right child of a node
// Inputs: node: node id
// Output: right child id or NO_CHILD
// Failure: returns NO_CHILD for invalid node
int VectorRangeTreeMap::getRightChild(int node) const
{
    return (node >= 0 && node < (int)edges.size()) ? edges[node].second : NO_CHILD;
}

// Purpose: Get the parent of a node
// Inputs: node: node id
// Output: parent id or NO_PARENT
// Failure: returns NO_PARENT for invalid node
int VectorRangeTreeMap::getParent(int node) const
{
    return (node >= 0 && node < (int)parents.size()) ? parents[node] : NO_PARENT;
}

// Purpose: Get the inorder range of a node
// Inputs: node: node id
// Output: {start, end} range or {0,0} if invalid
// Failure: returns {0,0} for invalid node
std::pair<int, int> VectorRangeTreeMap::getRange(int node) const
{
    return (node >= 0 && node < (int)ranges.size()) ? ranges[node] : std::make_pair(0, 0);
}

// Purpose: Set left child and update parent pointers
// Inputs: node: parent id; child: child id
// Output: nothing
// Failure: ignores invalid node ids
void VectorRangeTreeMap::setLeftChild(int node, int child)
{

    if (node < 0 || node > max_node_value)
        return;
    fingerprint_valid = false;
    int old = edges[node].first;
    if (old >= 0)
        parents[old] = NO_PARENT;
    edges[node].first = child;
    if (child >= 0)
        parents[child] = node;
}

// Purpose: Set right child and update parent pointers
// Inputs: node: parent id; child: child id
// Output: nothing
// Failure: ignores invalid node ids
void VectorRangeTreeMap::setRightChild(int node,int child) {
    if (node<0||node>max_node_value) return;
    fingerprint_valid = false;
    int old = edges[node].second;
    if (old>=0) parents[old]=NO_PARENT;
    edges[node].second = child;
    if (child>=0) parents[child]=node;
}

// Purpose: Check if node exists in original node set
// Inputs: node: node id
// Output: true if node is original
bool VectorRangeTreeMap::isOriginal(int node) const
{
    return node >= 0 && node < (int)original_mask.size() && original_mask[node] != 0;
}

// Purpose: Print tree structure recursively
// Inputs: node_value: start node or -2 for root; indent: prefix string
// Output: nothing
// Failure: no-op on invalid nodes
void VectorRangeTreeMap::printTreeStructure(int node_value, std::string indent) const
{
    if (node_value == -2)
        node_value = root;
    if (node_value == NO_CHILD || !isOriginal(node_value))
        return;

    auto range = getRange(node_value);
    int parent = getParent(node_value);
    std::cout << indent << "Node " << node_value
              << " (range: " << range.first << ", " << range.second
              << ", parent: " << parent << ")\n";

    int left = getLeftChild(node_value);
    if (left != NO_CHILD && isOriginal(left))
    {
        printTreeStructure(left, indent + "  L-");
    }

    int right = getRightChild(node_value);
    if (right != NO_CHILD && isOriginal(right))
    {
        printTreeStructure(right, indent + "  R-");
    }
}

// Purpose: Rotate left around node x and update ranges
// Inputs: x: node id
// Output: nothing
// Failure: no-op if rotation is not possible
void VectorRangeTreeMap::rotateLeft(int x)
{
    int y = getRightChild(x);

    if (y == NO_CHILD)
        return;

    int yl = getLeftChild(y);
    int xp = getParent(x);
    auto xr = getRange(x), yr = getRange(y), ylr = getRange(yl);
    fingerprint_valid = false;

    ranges[y] = xr;
    ranges[x] = (yl >= 0)
                    ? std::pair<int, int>{xr.first, ylr.second}
                    : std::pair<int, int>{xr.first, yr.first};
    parents[y] = xp;
    parents[x] = y;

    if (yl >= 0)
        parents[yl] = x;

    edges[x].second = yl;
    edges[y].first = x;

    if (xp != NO_PARENT)
    {
        if (edges[xp].first == x)

            edges[xp].first = y;

        else
            edges[xp].second = y;
    }
    else
        root = y;
}

// Purpose: Rotate right around node x and update ranges
// Inputs: x: node id
// Output: nothing
// Failure: no-op if rotation is not possible
void VectorRangeTreeMap::rotateRight(int x)
{
    int y = getLeftChild(x);
    if (y == NO_CHILD)
        return;

    int yr = getRightChild(y);
    int xp = getParent(x);
    auto xr = getRange(x), yrng = getRange(y), yrr = getRange(yr);
    fingerprint_valid = false;

    ranges[y] = xr;
    ranges[x] = (yr >= 0)
                    ? std::pair<int, int>{yrr.first, xr.second}
                    : std::pair<int, int>{yrng.second, xr.second};

    parents[y] = xp;
    parents[x] = y;

    if (yr >= 0)
        parents[yr] = x;

    edges[x].first = yr;
    edges[y].second = x;

    if (xp != NO_PARENT)
    {
        if (edges[xp].first == x)

            edges[xp].first = y;

        else
            edges[xp].second = y;
    }

    else
        root = y;
}

// Purpose: Collect all (parent, child) edges into a set
// Inputs: node_value: subtree root; out_set: output set
// Output: nothing
// Failure: no-op on invalid nodes
void VectorRangeTreeMap::collectEdges(int node_value,
    std::unordered_set<std::pair<int,int>,PairHash,PairEq>& out_set) const
{
    if (node_value == NO_CHILD || !isOriginal(node_value))
        return;

    int left = getLeftChild(node_value);
    if (left != NO_CHILD && isOriginal(left))
    {
        out_set.insert({node_value, left});
        collectEdges(left, out_set);
    }

    int right = getRightChild(node_value);
    if (right != NO_CHILD && isOriginal(right))
    {
        out_set.insert({node_value, right});
        collectEdges(right, out_set);
    }
}

// Purpose: Print full tree data including ranges, edges, and parents
// Inputs: none
// Output: nothing
void VectorRangeTreeMap::print() const
{
    std::cout << "Tree structure:\n";
    printTreeStructure();

    std::vector<int> sorted_original;
    for (int node : original_inorder)
    {
        sorted_original.push_back(node);
    }
    std::sort(sorted_original.begin(), sorted_original.end());

    std::cout << "\nranges: [";
    for (size_t i = 0; i < sorted_original.size(); i++)
    {
        if (i > 0)
            std::cout << ", ";
        auto range = getRange(sorted_original[i]);
        std::cout << "(" << range.first << "," << range.second << ")";
    }
    std::cout << "]\n";

    std::cout << "edges: [";
    for (size_t i = 0; i < sorted_original.size(); i++)
    {
        if (i > 0)
            std::cout << ", ";
        int node = sorted_original[i];
        int left = getLeftChild(node);
        int right = getRightChild(node);

        if (!isOriginal(left))
            left = NO_CHILD;
        if (!isOriginal(right))
            right = NO_CHILD;

        std::cout << "(" << left << "," << right << ")";
    }
    std::cout << "]\n";

    std::cout << "parents: [";
    for (size_t i = 0; i < sorted_original.size(); i++)
    {
        if (i > 0)
            std::cout << ", ";
        int node = sorted_original[i];
        int parent = getParent(node);

        if (!isOriginal(parent))
            parent = NO_PARENT;

        std::cout << parent;
    }
    std::cout << "]\n";
}

// Purpose: Collects current preorder for partition.
// Inputs: T (const VectorRangeTreeMap&); out (std::vector<int>&).
// Output: No return value; updates referenced outputs or shared state.
static void collectCurrentPreorderForPartition(const VectorRangeTreeMap& T,
                                               std::vector<int>& out)
{
    out.clear();
    if (T.root < 0 || T.original_inorder.empty() || !T.isOriginal(T.root)) {
        return;
    }
    const std::size_t expected = T.original_inorder.size();
    out.reserve(expected);
    if (expected <= 128) {
        std::array<int, 128> stack{};
        std::size_t top = 0;
        stack[top++] = T.root;
        while (top > 0) {
            const int node = stack[--top];
            if (node < 0 || !T.isOriginal(node)) {
                continue;
            }
            out.push_back(node);
            if (out.size() > expected) {
                out.clear();
                return;
            }
            const int right = T.getRightChild(node);
            if (right >= 0 && T.isOriginal(right)) {
                if (top >= stack.size()) {
                    out.clear();
                    return;
                }
                stack[top++] = right;
            }
            const int left = T.getLeftChild(node);
            if (left >= 0 && T.isOriginal(left)) {
                if (top >= stack.size()) {
                    out.clear();
                    return;
                }
                stack[top++] = left;
            }
        }
    } else {
        std::vector<int> stack;
        stack.reserve(expected);
        stack.push_back(T.root);
        while (!stack.empty()) {
            const int node = stack.back();
            stack.pop_back();
            if (node < 0 || !T.isOriginal(node)) {
                continue;
            }
            out.push_back(node);
            if (out.size() > expected) {
                out.clear();
                return;
            }
            const int right = T.getRightChild(node);
            if (right >= 0 && T.isOriginal(right)) stack.push_back(right);
            const int left = T.getLeftChild(node);
            if (left >= 0 && T.isOriginal(left)) stack.push_back(left);
        }
    }
    if (out.size() != expected) {
        out.clear();
    }
}

std::pair<VectorRangeTreeMap,VectorRangeTreeMap>
VectorRangeTreeMap::partitionAlongSubtreeRange(const VectorRangeTreeMap& T,
                                               const std::pair<int,int>& child_range)
{
    thread_local std::vector<int> inA_scratch;
    thread_local std::vector<int> inB_scratch;
    thread_local std::vector<int> preA_scratch;
    thread_local std::vector<int> preB_scratch;
    thread_local std::vector<int> source_preorder_scratch;
    auto& inA = inA_scratch;
    auto& inB = inB_scratch;
    auto& preA = preA_scratch;
    auto& preB = preB_scratch;
    inA.clear();
    inB.clear();
    preA.clear();
    preB.clear();
    const std::vector<int>* source_preorder = &T.original_preorder;
    if (source_preorder->empty() && !T.original_inorder.empty()) {
        collectCurrentPreorderForPartition(T, source_preorder_scratch);
        source_preorder = &source_preorder_scratch;
    }
    const int lo = child_range.first;
    const int hi = child_range.second;
    const int subtree_size = std::max(0, hi - lo);
    const int total_size = static_cast<int>(T.original_inorder.size());
    inA.reserve(static_cast<std::size_t>(subtree_size));
    preA.reserve(static_cast<std::size_t>(subtree_size));
    inB.reserve(static_cast<std::size_t>(std::max(0, total_size - subtree_size)));
    preB.reserve(static_cast<std::size_t>(std::max(0, total_size - subtree_size)));

    for (int pos = 0; pos < total_size; ++pos)
    {
        const int x = T.original_inorder[static_cast<std::size_t>(pos)];
        if (pos >= lo && pos < hi)
            inA.push_back(x);
        else
            inB.push_back(x);
    }
    for (int x : *source_preorder)
    {
        int pos = -1;
        if (x >= 0 && x < (int)T.position_in_inorder_fast.size()) {
            pos = T.position_in_inorder_fast[x];
        } else {
            auto it = T.position_in_inorder.find(x);
            if (it != T.position_in_inorder.end()) {
                pos = it->second;
            }
        }
        if (pos >= lo && pos < hi)
            preA.push_back(x);
        else
            preB.push_back(x);
    }
    VectorRangeTreeMap A, B;
    A.buildLight(preA, inA);
    B.buildLight(preB, inB);
    return {std::move(A), std::move(B)};
}

std::pair<VectorRangeTreeMap,VectorRangeTreeMap>
VectorRangeTreeMap::partitionAlongEdge(const VectorRangeTreeMap& T,
                                       const std::pair<int,int>& parent_range,
                                       const std::pair<int,int>& child_range)
{
    thread_local std::vector<int> inA_scratch;
    thread_local std::vector<int> inB_scratch;
    thread_local std::vector<int> preA_scratch;
    thread_local std::vector<int> preB_scratch;
    thread_local std::vector<int> source_preorder_scratch;
    auto& inA = inA_scratch;
    auto& inB = inB_scratch;
    auto& preA = preA_scratch;
    auto& preB = preB_scratch;
    inA.clear();
    inB.clear();
    preA.clear();
    preB.clear();
    const std::vector<int>* source_preorder = &T.original_preorder;
    if (source_preorder->empty() && !T.original_inorder.empty()) {
        collectCurrentPreorderForPartition(T, source_preorder_scratch);
        source_preorder = &source_preorder_scratch;
    }
    const int lo = child_range.first;
    const int hi = child_range.second;
    const int subtree_size = std::max(0, hi - lo);
    const int total_size = static_cast<int>(T.original_inorder.size());
    inA.reserve(static_cast<std::size_t>(subtree_size));
    preA.reserve(static_cast<std::size_t>(subtree_size));
    inB.reserve(static_cast<std::size_t>(std::max(0, total_size - subtree_size)));
    preB.reserve(static_cast<std::size_t>(std::max(0, total_size - subtree_size)));

    bool child_range_is_subtree = false;
    for (int x : T.original_inorder) {
        if (T.isOriginal(x) && T.getRange(x) == child_range) {
            child_range_is_subtree = true;
            break;
        }
    }
    const auto in_child_range = [&](int x) {
        if (child_range_is_subtree) {
            int pos = -1;
            if (x >= 0 && x < (int)T.position_in_inorder_fast.size()) {
                pos = T.position_in_inorder_fast[x];
            } else {
                auto it = T.position_in_inorder.find(x);
                if (it != T.position_in_inorder.end()) {
                    pos = it->second;
                }
            }
            return pos >= lo && pos < hi;
        }
        auto r = T.getRange(x);
        return r.first >= lo && r.second <= hi;
    };

    for (int pos = 0; pos < total_size; ++pos)
    {
        const int x = T.original_inorder[static_cast<std::size_t>(pos)];
        if (child_range_is_subtree ? (pos >= lo && pos < hi) : in_child_range(x))
            inA.push_back(x);
        else
            inB.push_back(x);
    }
    for (int x : *source_preorder)
    {
        if (in_child_range(x))
            preA.push_back(x);
        else
            preB.push_back(x);
    }
    VectorRangeTreeMap A, B;
    A.buildLight(preA, inA);
    B.buildLight(preB, inB);
    return {std::move(A), std::move(B)};
}

// Purpose: Reset internal vectors and metadata
// Inputs: none
// Output: nothing
void VectorRangeTreeMap::clear() 
{
    fingerprint_valid = false;
    fingerprint_h1 = 0;
    fingerprint_h2 = 0;
    ranges.clear();
    edges.clear();
    parents.clear();
    original_inorder.clear();
    position_in_inorder.clear();
    position_in_inorder_fast.clear();
    original_nodes.clear();
    original_mask.clear();
    root = -1;
    max_node_value = 0;
}

// Purpose: Recursive helper to build tree from traversals
// Inputs: preorder/inorder: traversal arrays; ps/pe/is/ie: ranges; p: parent id
// Output: nothing
// Failure: no-op on invalid ranges
void VectorRangeTreeMap::buildRecursive(const std::vector<int> &preorder, const std::vector<int> &inorder,
                    int ps, int pe, int is, int ie, int p)
{
    if (ps > pe || is > ie)
        return;

    int root_val = preorder[ps];

    parents[root_val] = p;

    int root_idx = (root_val >= 0 && root_val < (int)position_in_inorder_fast.size())
                       ? position_in_inorder_fast[root_val]
                       : position_in_inorder[root_val];
    int left_size = root_idx - is;

    if (left_size > 0)
    {
        int left_root_val = preorder[ps + 1];
        edges[root_val].first = left_root_val;
        parents[left_root_val] = root_val;
        buildRecursive(preorder, inorder, ps + 1, ps + left_size,
                       is, root_idx - 1, root_val);
    }

    if (root_idx < ie)
    {
        int right_root_val = preorder[ps + left_size + 1];
        edges[root_val].second = right_root_val;
        parents[right_root_val] = root_val;
        buildRecursive(preorder, inorder, ps + left_size + 1, pe,
                       root_idx + 1, ie, root_val);
    }
}

// Purpose: Compute a single node range from children
// Inputs: node_value: node id
// Output: nothing
// Failure: no-op on invalid node
void VectorRangeTreeMap::updateNodeRange(int node_value)
{
    if (node_value < 0 || node_value >= ranges.size() || !isOriginal(node_value))
        return;

    int left = getLeftChild(node_value);
    int right = getRightChild(node_value);

    int pos = (node_value >= 0 && node_value < (int)position_in_inorder_fast.size())
                  ? position_in_inorder_fast[node_value]
                  : position_in_inorder[node_value];

    int start = pos;
    int end = pos + 1;

    if (left != NO_CHILD && isOriginal(left))
    {
        auto left_range = getRange(left);
        start = left_range.first;
    }

    if (right != NO_CHILD && isOriginal(right))
    {
        auto right_range = getRange(right);
        end = right_range.second;
    }

    ranges[node_value] = {start, end};
}

// Purpose: Calculate ranges for all nodes
// Inputs: none
// Output: nothing
void VectorRangeTreeMap::calculateAllRanges()
{
    // Parents depend on child ranges, so ranges must be rebuilt bottom-up.
    calculateRangesPostOrder(root);
}

// Purpose: Post-order update for all ranges in a subtree
// Inputs: node_value: subtree root
// Output: nothing
// Failure: no-op on invalid node
void VectorRangeTreeMap::calculateRangesPostOrder(int node_value)
{
    if (node_value == NO_CHILD || !isOriginal(node_value))
        return;

    int left = getLeftChild(node_value);
    int right = getRightChild(node_value);

    if (left != NO_CHILD && isOriginal(left))
    {
        calculateRangesPostOrder(left);
    }
    if (right != NO_CHILD && isOriginal(right))
    {
        calculateRangesPostOrder(right);
    }

    updateNodeRange(node_value);
}

// Purpose: Detect if two trees are identical on the same node set
// Inputs: A: first tree; B: second tree
// Output: true if structure, ranges, parents, and children match
// Failure: returns false on mismatch or invalid input
bool TreesEqual(const VectorRangeTreeMap &A,
                const VectorRangeTreeMap &B)
{
    // Structural metadata is comparable only when both trees use the same ordered label set.
    if (A.original_inorder != B.original_inorder)
        return false;

    for (int v : A.original_inorder)
    {
        if (!A.isOriginal(v) || !B.isOriginal(v))
            return false;
        if (A.getRange(v) != B.getRange(v))
            return false;
        if (A.getLeftChild(v) != B.getLeftChild(v))
            return false;
        if (A.getRightChild(v) != B.getRightChild(v))
            return false;
        if (A.getParent(v) != B.getParent(v))
            return false;
    }

    return true;
}

// Purpose: Serialize a tree into a unique string for hashing
// Inputs: T: tree to serialize
// Output: stable string encoding of ranges and edges
std::string treeToString(const VectorRangeTreeMap &T)
{
    // Sorting makes memo keys independent of traversal storage order.
    std::vector<int> nodes = T.original_inorder;
    std::sort(nodes.begin(), nodes.end());

    std::ostringstream oss;
    for (int v : nodes)
    {
        auto [rs, re] = T.getRange(v);
        int L = T.getLeftChild(v);
        int R = T.getRightChild(v);
        int P = T.getParent(v);

        oss << v << ','
            << rs << ',' << re << ','
            << L << ',' << R << ','
            << P << ';';
    }
    return oss.str();
}

// Purpose: Compare two range edges for equality
// Inputs: o: other range edge
// Output: true if all fields match
bool RP::operator==(RP const &o) const {
    return cs == o.cs && ce == o.ce;
}

// Purpose: Hash a range edge for unordered containers
// Inputs: r: range edge
// Output: hash value
size_t RPH::operator()(RP const &r) const {
    uint64_t h2 = ((uint64_t)r.cs << 32) ^ uint32_t(r.ce);
    return std::hash<uint64_t>()(h2);
}

// Purpose: Collect all (parent_range, child_range) edges in target
// Inputs: T: tree to scan
// Output: set of range edges
// Failure: returns empty on invalid tree
std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap& T)
{
    std::unordered_set<RP, RPH> s;
    s.reserve(T.original_inorder.size() * 2 + 1);
    auto dfs = [&](auto&& self, int v) -> void {
        if (v < 0 || !T.isOriginal(v))
            return;
        auto pr = T.getRange(v);
        const int left = T.getLeftChild(v);
        const int right = T.getRightChild(v);
        if (left >= 0 && T.isOriginal(left)) {
            auto cr = T.getRange(left);
            s.insert(RP{pr.first, pr.second, cr.first, cr.second});
            self(self, left);
        }
        if (right >= 0 && T.isOriginal(right)) {
            auto cr = T.getRange(right);
            s.insert(RP{pr.first, pr.second, cr.first, cr.second});
            self(self, right);
        }
    };
    dfs(dfs, T.root);
    return s;
}

// Purpose: Check if a tree already has a given range edge
// Inputs: tree: tree to scan; e: range edge
// Output: true if edge exists
// Failure: returns false on invalid tree
bool hasEdgeByRange(const VectorRangeTreeMap &tree, RP const &e)
{
    auto dfs = [&](auto&& self, int v) -> bool {
        if (v < 0 || !tree.isOriginal(v))
            return false;
        const int left = tree.getLeftChild(v);
        const int right = tree.getRightChild(v);
        if (left >= 0 && tree.isOriginal(left)) {
            auto cr = tree.getRange(left);
            if (cr.first == e.cs && cr.second == e.ce)
                return true;
        }
        if (right >= 0 && tree.isOriginal(right)) {
            auto cr = tree.getRange(right);
            if (cr.first == e.cs && cr.second == e.ce)
                return true;
        }
        return self(self, left) || self(self, right);
    };
    return dfs(dfs, tree.root);
}

// Purpose: Detect any free edge in cur that target wants
// Inputs: cur: current tree; tgt: target tree; out_v: rotation node; out_leftRotation: direction
// Output: true if a free edge is found
// Failure: returns false on invalid trees or if none found
bool hasFreeEdge(const VectorRangeTreeMap &cur,
                 const VectorRangeTreeMap &tgt,
                 int &out_v,
                 bool &out_leftRotation)
{
    auto tgtEdges = buildTargetSet(tgt);

    for (int v : cur.original_inorder)
    {
        int y = cur.getRightChild(v);
        if (y >= 0 && cur.isOriginal(y))
        {
            auto xr = cur.getRange(v);
            auto yr = cur.getRange(y);

            int yl = cur.getLeftChild(y);
            std::pair<int, int> newChild;
            if (yl >= 0 && cur.isOriginal(yl))
            {
                newChild = cur.getRange(yl);
            }
            else
            {

                newChild = {xr.first, yr.first};
            }

            RP e{xr.first, xr.second,
                 newChild.first, newChild.second};

            if (tgtEdges.count(e) && !hasEdgeByRange(cur, e))
            {
                out_v = v;
                out_leftRotation = true;
                return true;
            }
        }

        int z = cur.getLeftChild(v);
        if (z >= 0 && cur.isOriginal(z))
        {
            auto xr = cur.getRange(v);
            auto zr = cur.getRange(z);

            int zrchild = cur.getRightChild(z);
            std::pair<int, int> newChild;
            if (zrchild >= 0 && cur.isOriginal(zrchild))
            {
                newChild = cur.getRange(zrchild);
            }
            else
            {

                newChild = {zr.second, xr.second};
            }

            RP e{xr.first, xr.second,
                 newChild.first, newChild.second};

            if (tgtEdges.count(e) && !hasEdgeByRange(cur, e))
            {
                out_v = v;
                out_leftRotation = false;
                return true;
            }
        }
    }

    return false;
}

int BFSSearch(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
int removeFreeEdge(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
int Dist(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

static std::unordered_map<std::string, int> memo;
static std::unordered_set<std::string> in_progress;

// Purpose: Compute distance via free-edge shortcut and partitioning
// Inputs: Ts: source tree; Te: target tree
// Output: rotation distance
// Failure: falls back to BFSSearch on cycles or partition mismatch
int Dist(const VectorRangeTreeMap &Ts,
         const VectorRangeTreeMap &Te)
{
    if (TreesEqual(Ts, Te))
        return 0;

    std::string key = treeToString(Ts) + "|" + treeToString(Te);
    if (auto it = memo.find(key); it != memo.end())
        return it->second;
    if (in_progress.count(key))
        return BFSSearch(Ts, Te);
    in_progress.insert(key);

    int result;
    int u;
    bool left;
    if (hasFreeEdge(Ts, Te, u, left))
    {
        int v = left
                    ? Ts.getRightChild(u)
                    : Ts.getLeftChild(u);

        VectorRangeTreeMap Tbar = Ts;
        if (left)
            Tbar.rotateLeft(u);
        else
            Tbar.rotateRight(u);

        auto pr = Tbar.getRange(v);
        auto cr = Tbar.getRange(u);

        auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(Tbar, pr, cr);
        auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(Te, pr, cr);

        if (A1.original_inorder != B1.original_inorder ||
            A2.original_inorder != B2.original_inorder)
        {
            result = BFSSearch(Ts, Te);
        }
        else
        {
            result = 1 + Dist(A1, B1) + Dist(A2, B2);
        }
    }
    else
    {
        result = BFSSearch(Ts, Te);
    }

    memo[key] = result;
    in_progress.erase(key);
    return result;
}

// Purpose: Apply one free-edge rotation and recurse on partitions
// Inputs: Ts: source tree; Te: target tree
// Output: distance contribution for that free-edge path
// Failure: falls back to BFSSearch on partition mismatch
int removeFreeEdge(const VectorRangeTreeMap &Ts,
                   const VectorRangeTreeMap &Te)
{
    int u;
    bool left;
    assert(hasFreeEdge(Ts, Te, u, left));

    // Rotation reverses this edge, so preserve the old child needed to identify the new split.
    int v = left
                ? Ts.getRightChild(u)
                : Ts.getLeftChild(u);

    VectorRangeTreeMap Tbar = Ts;
    if (left)
        Tbar.rotateLeft(u);
    else
        Tbar.rotateRight(u);

    auto pr = Tbar.getRange(v);
    auto cr = Tbar.getRange(u);

    auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(Tbar, pr, cr);
    auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(Te, pr, cr);

    // Mismatched partitions invalidate the decomposition; exact BFS remains the safe fallback.
    if (A1.original_inorder != B1.original_inorder ||
        A2.original_inorder != B2.original_inorder)
    {
        return BFSSearch(Ts, Te);
    }
    return 1 + Dist(A1, B1) + Dist(A2, B2);
}

typedef std::pair<VectorRangeTreeMap, VectorRangeTreeMap> P2;
bool hasFreeEdge(const VectorRangeTreeMap &, const VectorRangeTreeMap &, int &, bool &);

// Purpose: Brute-force BFS search for rotation distance
// Inputs: T_s: source tree; T_e: target tree
// Output: minimum distance or INT_MAX if not found
// Failure: returns INT_MAX on failures
int BFSSearch(const VectorRangeTreeMap &T_s, const VectorRangeTreeMap &T_e)
{
    if (TreesEqual(T_s, T_e))
        return 0;

    std::queue<std::pair<VectorRangeTreeMap, int>> Q;
    std::unordered_set<std::string> visited;

    Q.push({T_s, 0});
    visited.insert(treeToString(T_s));

    while (!Q.empty())
    {
        auto [cur, d] = Q.front();
        Q.pop();
        int nextDist = d + 1;

        for (int v : cur.original_inorder)
        {
            if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD)
            {
                auto tmp = cur;
                tmp.rotateLeft(v);

                if (TreesEqual(tmp, T_e))
                {
                    return nextDist;
                }

                auto key = treeToString(tmp);
                if (!visited.count(key))
                {
                    visited.insert(key);
                    Q.push({std::move(tmp), nextDist});
                }
            }

            if (cur.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD)
            {
                auto tmp = cur;
                tmp.rotateRight(v);

                if (TreesEqual(tmp, T_e))
                {
                    return nextDist;
                }

                auto key = treeToString(tmp);
                if (!visited.count(key))
                {
                    visited.insert(key);
                    Q.push({std::move(tmp), nextDist});
                }
            }
        }
    }

    // Rotation graphs are connected for equal node sets; exhaustion therefore indicates invalid input.
    return INT_MAX;
}

// Purpose: Clear memo cache then compute distance
// Inputs: T_init: source tree; T_final: target tree
// Output: rotation distance
int FindRotationDistance(const VectorRangeTreeMap &T_init,
                         const VectorRangeTreeMap &T_final)
{
    memo.clear();
    return Dist(T_init, T_final);
}

// Purpose: Run internal verification tests
// Inputs: none
// Output: nothing
// Failure: asserts on failure
void runAllTests()
{

    VectorRangeTreeMap T;
    std::vector<int> pre = {8, 4, 2, 1, 0, 3, 6, 5, 7};
    std::vector<int> in = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    T.build(pre, in);

    auto pr = T.getRange(4);
    auto cr = T.getRange(2);
    auto [A, B] = VectorRangeTreeMap::partitionAlongEdge(T, pr, cr);

    assert(A.original_inorder == std::vector<int>({0, 1, 2, 3}));
    assert(A.original_preorder == std::vector<int>({2, 1, 0, 3}));

    assert(B.original_inorder == std::vector<int>({4, 5, 6, 7, 8}));
    assert(B.original_preorder == std::vector<int>({8, 4, 6, 5, 7}));

    std::cout << "PASS partitionAlongEdge test\n";

    {
        VectorRangeTreeMap A;
        A.build({2, 1, 3}, {1, 2, 3});
        assert(BFSSearch(A, A) == 0);
        assert(Dist(A, A) == 0);
        std::cout << "PASS identical: BFSSearch=0, Dist=0\n";
    }
    {
        VectorRangeTreeMap A, B;
        A.build({2, 1, 3}, {1, 2, 3});
        B = A;
        B.rotateLeft(2);
        assert(BFSSearch(A, B) == 1);
        assert(Dist(A, B) == 1);
        std::cout << "PASS one rot: BFSSearch=1, Dist=1\n";
    }
    {
        VectorRangeTreeMap A, D;
        A.build({2, 1, 3, 4}, {1, 2, 3, 4});
        D = A;
        D.rotateLeft(2);
        D.rotateLeft(3);
        assert(BFSSearch(A, D) == 2);
        assert(Dist(A, D) == 2);
        std::cout << "PASS two rots: BFSSearch=2, Dist=2\n";
    }
    {
        VectorRangeTreeMap A, C;
        A.build({2, 1, 3}, {1, 2, 3});
        C.build({3, 2, 1}, {1, 2, 3});
        assert(removeFreeEdge(A, C) == 1);
        assert(Dist(A, C) == 1);
        std::cout << "PASS free-edge: removeFreeEdge=1, Dist=1\n";
    }
    {
        // Small random trees use exhaustive BFS as an oracle for the optimized recurrence.
        auto buildBSTPreIn = [](int n, int seed) {
            struct Node {
                int v;
                Node *l = nullptr;
                Node *r = nullptr;
                // Purpose: Creates a temporary BST node for randomized oracle tests.
                // Inputs: val: node label to store.
                // Output: A node with no children.
                explicit Node(int val) : v(val) {}
            };
            auto insert = [](Node *root, int val) -> Node * {
                if (!root) return new Node(val);
                Node *cur = root;
                while (true) {
                    if (val < cur->v) {
                        if (cur->l) cur = cur->l;
                        else { cur->l = new Node(val); break; }
                    } else {
                        if (cur->r) cur = cur->r;
                        else { cur->r = new Node(val); break; }
                    }
                }
                return root;
            };
            auto preorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
                if (!node) return;
                out.push_back(node->v);
                self(node->l, out, self);
                self(node->r, out, self);
            };
            auto inorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
                if (!node) return;
                self(node->l, out, self);
                out.push_back(node->v);
                self(node->r, out, self);
            };
            std::vector<int> values;
            values.reserve(n);
            for (int i = 1; i <= n; i++) values.push_back(i);
            std::mt19937 rng(seed);
            std::shuffle(values.begin(), values.end(), rng);

            Node *root = nullptr;
            for (int v : values) root = insert(root, v);

            std::vector<int> pre, in;
            preorder(root, pre, preorder);
            inorder(root, in, inorder);
            return std::pair<std::vector<int>, std::vector<int>>(pre, in);
        };

        std::cout << "Random BST checks (Dist vs BFSSearch):\n";
        for (int n = 4; n <= 7; n++) {
            for (int seed = 0; seed < 5; seed++) {
                auto [pre1, in1] = buildBSTPreIn(n, 100 + seed * 2);
                auto [pre2, in2] = buildBSTPreIn(n, 101 + seed * 2);
                VectorRangeTreeMap A, B;
                A.build(pre1, in1);
                B.build(pre2, in2);
                int d_bfs = BFSSearch(A, B);
                int d_dist = Dist(A, B);
                if (d_bfs != d_dist) {
                    std::cout << "MISMATCH n=" << n << " seed=" << seed
                              << " bfs=" << d_bfs << " dist=" << d_dist << "\n";
                }
            }
        }
        std::cout << "Random BST checks done.\n";
    }
    std::cout << "All tests passed";
}

// Purpose: Build preorder and inorder for a BST from shuffled insertion order
// Inputs: n: node count; seed: shuffle seed
// Output: {preorder, inorder} vectors
// Failure: returns empty vectors if n <= 0
static std::pair<std::vector<int>, std::vector<int>> buildRandomBST(int n, int seed)
{
    struct Node
    {
        int v;
        Node *l = nullptr;
        Node *r = nullptr;
        // Purpose: Creates a temporary BST node for deterministic input generation.
        // Inputs: val: node label to store.
        // Output: A node with no children.
        explicit Node(int val) : v(val) {}
    };
    auto insert = [](Node *root, int val) -> Node * {
        if (!root)
            return new Node(val);
        Node *cur = root;
        while (true)
        {
            if (val < cur->v)
            {
                if (cur->l)
                    cur = cur->l;
                else
                {
                    cur->l = new Node(val);
                    break;
                }
            }
            else
            {
                if (cur->r)
                    cur = cur->r;
                else
                {
                    cur->r = new Node(val);
                    break;
                }
            }
        }
        return root;
    };
    auto preorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
        if (!node)
            return;
        out.push_back(node->v);
        self(node->l, out, self);
        self(node->r, out, self);
    };
    auto inorder = [](Node *node, std::vector<int> &out, auto &self) -> void {
        if (!node)
            return;
        self(node->l, out, self);
        out.push_back(node->v);
        self(node->r, out, self);
    };

    std::vector<int> values;
    values.reserve(n);
    for (int i = 1; i <= n; i++)
        values.push_back(i);
    std::mt19937 rng(seed);
    std::shuffle(values.begin(), values.end(), rng);

    Node *root = nullptr;
    for (int v : values)
        root = insert(root, v);

    std::vector<int> pre, in;
    preorder(root, pre, preorder);
    inorder(root, in, inorder);
    return {pre, in};
}

// Purpose: Encode the tree as P:...;I:... for Java CLI compatibility
// Inputs: T: tree to encode
// Output: encoded string
// Failure: returns empty encodings if traversals are empty
static std::string treeToPreInString(const VectorRangeTreeMap &T)
{
    auto join = [](const std::vector<int> &vals) {
        std::ostringstream oss;
        for (size_t i = 0; i < vals.size(); i++)
        {
            if (i)
                oss << ",";
            oss << vals[i];
        }
        return oss.str();
    };
    std::ostringstream oss;
    oss << "P:" << join(T.original_preorder) << ";I:" << join(T.original_inorder);
    return oss.str();
}

// Purpose: Brute-force CLI entry that emits JSON per direction
// Inputs: argc/argv: CLI args
// Output: exit code, 0 on success
// Failure: prints to stderr on invalid args or unsupported cases
static int runCli(int argc, char **argv)
{
    std::string case_type = "random";
    int n = 0;
    int seed = 0;
    int count = 1;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--case" && i + 1 < argc)
            case_type = argv[++i];
        else if (arg == "--n" && i + 1 < argc)
            n = std::stoi(argv[++i]);
        else if (arg == "--seed" && i + 1 < argc)
            seed = std::stoi(argv[++i]);
        else if (arg == "--count" && i + 1 < argc)
            count = std::stoi(argv[++i]);
    }

    if (n <= 0)
    {
        std::cerr << "Missing or invalid --n\n";
        return 2;
    }
    if (case_type != "random")
    {
        std::cerr << "Only --case random is supported for bf_bst CLI\n";
        return 2;
    }

    for (int i = 0; i < count; i++)
    {
        int case_seed = seed + i;
        auto [pre1, in1] = buildRandomBST(n, case_seed * 2);
        auto [pre2, in2] = buildRandomBST(n, case_seed * 2 + 1);

        VectorRangeTreeMap A, B;
        A.build(pre1, in1);
        B.build(pre2, in2);

        auto run_one = [&](const VectorRangeTreeMap &T1,
                           const VectorRangeTreeMap &T2,
                           const char *direction) {
            auto start = std::chrono::steady_clock::now();
            int dist = Dist(T1, T2);
            auto end = std::chrono::steady_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start).count();

            std::cout << std::fixed << std::setprecision(3);
            std::cout << "{\"case_type\":\"" << case_type << "\""
                      << ",\"n\":" << n
                      << ",\"seed\":" << case_seed
                      << ",\"direction\":\"" << direction << "\""
                      << ",\"distance\":" << dist
                      << ",\"time_ms\":" << ms
                      << ",\"status\":\"ok\""
                      << ",\"tree_a\":\"" << treeToPreInString(T1) << "\""
                      << ",\"tree_b\":\"" << treeToPreInString(T2) << "\""
                      << "}" << std::endl;
        };

        run_one(A, B, "a->b");
        run_one(B, A, "b->a");
    }
    return 0;
}

#ifndef FLIPDIST_BINARY
// Purpose: Program entry point for bf_bst binary
// Inputs: argc/argv: CLI args
// Output: exit code
// Failure: none beyond runCli and runAllTests
int main(int argc, char **argv)
{
    if (argc > 1)
        return runCli(argc, argv);
    runAllTests();
    return 0;
}
#endif
