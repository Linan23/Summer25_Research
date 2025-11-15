// Implements the core mechanics of `VectorRangeTreeMap`: building from
// traversals, maintaining parent/child metadata, verifying structure, and
// computing range information. Rotation helpers live in tree_rotate.cpp and
// reuse the primitives defined here.

#include "rotation_tree.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <unordered_set>

// Default-initialised tree with no nodes; Build will populate structure.
VectorRangeTreeMap::VectorRangeTreeMap() : root(-1), max_node_value(0) {}

#ifndef NDEBUG
// Performs an expensive structural audit ensuring parent/child/range metadata
// is consistent; invoked after mutations when assertions are enabled.
void VectorRangeTreeMap::verify() const {
    for (int v : original_nodes) {
        assert(v >= 0 && v <= max_node_value);
        assert(position_in_inorder.count(v) && "missing inorder position");
    }
    for (int v : original_nodes) {
        int L = getLeftChild(v), R = getRightChild(v);
        if (L != NO_CHILD) { assert(isOriginal(L)); assert(getParent(L) == v); }
        if (R != NO_CHILD) { assert(isOriginal(R)); assert(getParent(R) == v); }
        int p = getParent(v);
        if (p != NO_PARENT) {
            assert(isOriginal(p));
            int PL = getLeftChild(p), PR = getRightChild(p);
            assert(PL == v || PR == v);
        }
    }
    if (root != NO_CHILD) assert(getParent(root) == NO_PARENT);

    std::unordered_set<int> seen;
    std::function<void(int)> dfs = [&](int x){
        if (x == NO_CHILD || !isOriginal(x) || seen.count(x)) return;
        seen.insert(x);
        dfs(getLeftChild(x));
        dfs(getRightChild(x));
    };
    dfs(root);
    assert(seen == original_nodes && "disconnected or stray nodes");

    for (int v : original_nodes) {
        auto r = getRange(v);
        assert(r.first < r.second);
        int pos = position_in_inorder.at(v);
        assert(r.first <= pos && pos < r.second);
    }
}
#endif

/*
void VectorRangeTreeMap::build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
    clear();
    if (preorder.empty()) return;

    original_preorder = preorder;
    original_inorder  = inorder;
    for (int val : preorder) original_nodes.insert(val);

    max_node_value = *std::max_element(preorder.begin(), preorder.end());
    max_node_value = std::max(max_node_value, *std::max_element(inorder.begin(), inorder.end()));
    ranges .assign(max_node_value+1, {0,0});
    edges  .assign(max_node_value+1, {NO_CHILD, NO_CHILD});
    parents.assign(max_node_value+1, NO_PARENT);

    for (size_t i = 0; i < inorder.size(); ++i)
        position_in_inorder[inorder[i]] = (int)i;

    root = preorder[0];
    buildRecursive(preorder, inorder, 0, (int)preorder.size()-1, 0, (int)inorder.size()-1, NO_PARENT);
    calculateAllRanges();
    verify();
}
*/

// Rebuilds the tree from preorder and inorder traversals; populates metadata
// (ranges, parents, edge tables) and verifies structural consistency.
void VectorRangeTreeMap::build(const std::vector<int>& preorder,
                               const std::vector<int>& inorder) {
    clear();
    if (preorder.empty()) return;
    if (preorder.size() != inorder.size())
        throw std::logic_error("build: preorder and inorder sizes differ");

    original_preorder = preorder;
    original_inorder  = inorder;

    original_nodes.clear();
    for (int v : inorder) original_nodes.insert(v);

    // sanity: preorder must be a permutation of inorder
#ifndef NDEBUG
    for (int v : preorder) {
        if (!original_nodes.count(v))
            throw std::logic_error("build: preorder has value not in inorder");
    }
#endif

    max_node_value = std::max(
        *std::max_element(preorder.begin(), preorder.end()),
        *std::max_element(inorder.begin(),  inorder.end())
    );
    ranges .assign(max_node_value + 1, {0,0});
    edges  .assign(max_node_value + 1, {NO_CHILD, NO_CHILD});
    parents.assign(max_node_value + 1, NO_PARENT);

    // refresh this every build
    position_in_inorder.clear();
    position_in_inorder.reserve(inorder.size());
    for (int i = 0; i < (int)inorder.size(); ++i)
        position_in_inorder[inorder[i]] = i;

    root = preorder[0];
    buildRecursive(preorder, inorder,
                   /*ps=*/0, /*pe=*/(int)preorder.size()-1,
                   /*is=*/0, /*ie=*/(int)inorder.size()-1,
                   /*p=*/NO_PARENT);

    calculateAllRanges();
    verify();
}

// Accessors 
int VectorRangeTreeMap::getLeftChild(int node)  const { return (node>=0 && node<(int)edges.size()) ? edges[node].first  : NO_CHILD; }
int VectorRangeTreeMap::getRightChild(int node) const { return (node>=0 && node<(int)edges.size()) ? edges[node].second : NO_CHILD; }
int VectorRangeTreeMap::getParent(int node)     const { return (node>=0 && node<(int)parents.size()) ? parents[node]    : NO_PARENT; }
std::pair<int,int> VectorRangeTreeMap::getRange(int node) const {
    return (node>=0 && node<(int)ranges.size()) ? ranges[node] : std::make_pair(0,0);
}
bool VectorRangeTreeMap::isOriginal(int node) const { return original_nodes.count(node) > 0; }

//  Mutators
void VectorRangeTreeMap::setLeftChild(int node, int child) {
    if (node < 0 || node > max_node_value) return;
    int old = edges[node].first;
    if (old >= 0) parents[old] = NO_PARENT;
    edges[node].first = child;
    if (child >= 0) parents[child] = node;
}
void VectorRangeTreeMap::setRightChild(int node, int child) {
    if (node < 0 || node > max_node_value) return;
    int old = edges[node].second;
    if (old >= 0) parents[old] = NO_PARENT;
    edges[node].second = child;
    if (child >= 0) parents[child] = node;
}

// internal
// Resets the tree to an empty state; used prior to rebuilding from traversals.
void VectorRangeTreeMap::clear() {
    ranges.clear(); edges.clear(); parents.clear();
    original_inorder.clear(); original_preorder.clear();
    position_in_inorder.clear(); original_nodes.clear();
    root = -1; max_node_value = 0;
}

/*
void VectorRangeTreeMap::buildRecursive(const std::vector<int>& preorder, const std::vector<int>& inorder,
                                        int ps, int pe, int is, int ie, int p) {
    if (ps > pe || is > ie) return;
    int root_val = preorder[ps];
    parents[root_val] = p;
    int root_idx = position_in_inorder[root_val];
    int left_size = root_idx - is;

    if (left_size > 0) {
        int left_root_val = preorder[ps+1];
        setLeftChild(root_val, left_root_val);
        buildRecursive(preorder, inorder, ps+1, ps+left_size, is, root_idx-1, root_val);
    }
    if (root_idx < ie) {
        int right_root_val = preorder[ps + left_size + 1];
        setRightChild(root_val, right_root_val);
        buildRecursive(preorder, inorder, ps+left_size+1, pe, root_idx+1, ie, root_val);
    }
}
*/

// Recursively materialises the subtree for the given preorder/inorder slices.
void VectorRangeTreeMap::buildRecursive(const std::vector<int>& preorder,
                                        const std::vector<int>& inorder,
                                        int ps, int pe, int is, int ie, int p) {
    if (ps > pe || is > ie) return;

    // Slice bounds must be valid and same length
    if (ps < 0 || pe >= (int)preorder.size() || is < 0 || ie >= (int)inorder.size())
        throw std::out_of_range("buildRecursive: slice out of bounds");

    int pre_len = pe - ps + 1;
    int in_len  = ie - is + 1;
    if (pre_len != in_len)
        throw std::logic_error("buildRecursive: preorder/inorder slice length mismatch");

    // Root and its inorder position
    const int root_val = preorder[ps];
    auto it = position_in_inorder.find(root_val);
    if (it == position_in_inorder.end())
        throw std::logic_error("buildRecursive: root not found in inorder map");

    const int root_idx = it->second;
    if (root_idx < is || root_idx > ie)
        throw std::logic_error("buildRecursive: root index outside inorder window");

    parents[root_val] = p;

    const int left_size  = root_idx - is;      // #nodes in left subtree
    const int right_size = ie - root_idx;      // #nodes in right subtree

    // Left subtree
    if (left_size > 0) {
        // Left preorder window is [ps+1, ps+left_size]
        const int lps = ps + 1;
        const int lpe = ps + left_size;
        if (lps > lpe || lpe > pe)
            throw std::logic_error("buildRecursive: left subtree preorder window invalid");

        const int left_root_val = preorder[lps];
        setLeftChild(root_val, left_root_val);

        buildRecursive(preorder, inorder,
                       /*ps=*/lps, /*pe=*/lpe,
                       /*is=*/is,  /*ie=*/root_idx - 1,
                       /*p=*/root_val);
    }

    // Right subtree
    if (right_size > 0) {
        // Right preorder window is [ps+left_size+1, pe]
        const int rps = ps + left_size + 1;
        const int rpe = pe;
        if (rps > rpe)
            throw std::logic_error("buildRecursive: right subtree preorder window invalid");

        const int right_root_val = preorder[rps];
        setRightChild(root_val, right_root_val);

        buildRecursive(preorder, inorder,
                       /*ps=*/rps,         /*pe=*/rpe,
                       /*is=*/root_idx+1,  /*ie=*/ie,
                       /*p=*/root_val);
    }
}

// Recomputes the inorder span for node v based on its current children.
void VectorRangeTreeMap::updateNodeRange(int v) {
    if (v < 0 || v >= (int)ranges.size() || !isOriginal(v)) return;
    int L = getLeftChild(v), R = getRightChild(v);
    int pos = position_in_inorder[v];
    int start = pos, end = pos+1;
    if (L != NO_CHILD && isOriginal(L)) start = getRange(L).first;
    if (R != NO_CHILD && isOriginal(R)) end   = getRange(R).second;
    ranges[v] = {start, end};
}

// Propagates range updates upward until the values stabilise.
void VectorRangeTreeMap::recomputeUpwardsFrom(int start) {
    int v = start;
    while (v != NO_PARENT) {
        auto old = ranges[v];
        updateNodeRange(v);
        if (ranges[v] == old) break;
        v = getParent(v);
    }
}

// Post-order traversal that seeds ranges for every node.
void VectorRangeTreeMap::calculateRangesPostOrder(int v) {
    if (v == NO_CHILD || !isOriginal(v)) return;
    int L = getLeftChild(v), R = getRightChild(v);
    if (L != NO_CHILD && isOriginal(L)) calculateRangesPostOrder(L);
    if (R != NO_CHILD && isOriginal(R)) calculateRangesPostOrder(R);
    updateNodeRange(v);
}
void VectorRangeTreeMap::calculateAllRanges() { calculateRangesPostOrder(root); }

void VectorRangeTreeMap::collectEdges(
    int node,
    std::unordered_set<std::pair<int,int>,PairHash,PairEq>& out) const
{
    if (node == NO_CHILD || !isOriginal(node)) return;
    int L = getLeftChild(node);
    if (L != NO_CHILD && isOriginal(L)) { out.insert({node,L}); collectEdges(L, out); }
    int R = getRightChild(node);
    if (R != NO_CHILD && isOriginal(R)) { out.insert({node,R}); collectEdges(R, out); }
}

std::pair<int,int> VectorRangeTreeMap::diagonalEndpoints(int node) const {
    if (node < 0 || node >= (int)ranges.size()) return {-1,-1};
    if (!isOriginal(node)) return {-1,-1};
    return ranges[node];
}

// Splits the tree along the parent->child edge defined by the ranges, returning
// the induced subtrees on either side of the cut.
std::pair<VectorRangeTreeMap, VectorRangeTreeMap>
VectorRangeTreeMap::partitionAlongEdge(const VectorRangeTreeMap& T,
    const std::pair<int,int>& /*parent_range*/,
    const std::pair<int,int>& child_range)
{
    std::vector<int> inA, inB, preA, preB;
    for (int x : T.original_inorder) {
        auto r = T.getRange(x);
        if (r.first >= child_range.first && r.second <= child_range.second) inA.push_back(x);
        else inB.push_back(x);
    }
    for (int x : T.original_preorder) {
        auto r = T.getRange(x);
        if (r.first >= child_range.first && r.second <= child_range.second) preA.push_back(x);
        else preB.push_back(x);
    }
    VectorRangeTreeMap A, B; A.build(preA, inA); B.build(preB, inB);
    return {A,B};
}
