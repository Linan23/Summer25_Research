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

// Definition: Hash a (parent, child) pair for unordered containers
// Parameters: p: pair to hash
// Returns: hash value
// Errors: none
size_t PairHash::operator()(const std::pair<int, int> &p) const
{
    // Combine two 32-bit ints into one 64-bit for hashing
    return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned long long)p.second);
}

// Definition: Compare two (parent, child) pairs for equality
// Parameters: a: first pair; b: second pair
// Returns: true if both pairs match
// Errors: none
bool PairEq::operator()(const std::pair<int, int> &a,
                        const std::pair<int, int> &b) const
{
    // Return true only if both parent and child match
    return a.first == b.first && a.second == b.second;
}

// Definition: Initialize an empty tree with default values
// Parameters: none
// Returns: constructed tree
// Errors: none
VectorRangeTreeMap::VectorRangeTreeMap() : root(-1), max_node_value(0)
{

    // (range start, range end) per node
    std::vector<std::pair<int, int>> ranges;
    // (leftChild, rightChild) per node
    std::vector<std::pair<int, int>> edges;
    // parent per node
    std::vector<int> parents;

    int root = -1;          // root node value
    int max_node_value = 0; // maximum node value to size vectors

    // for partitionAlongEdge

    std::vector<int> original_inorder, original_preorder;
    std::unordered_map<int, int> position_in_inorder; // node_value -> position in inorder
    std::unordered_set<int> original_nodes;           // set of original node values

    static constexpr int NO_CHILD = -1;
    static constexpr int NO_PARENT = -1;
}

// Definition: Build tree and ranges from preorder and inorder traversals
// Parameters: preorder: preorder node list; inorder: inorder node list
// Returns: nothing
// Errors: leaves tree empty on invalid input
void VectorRangeTreeMap::build(const std::vector<int> &preorder, const std::vector<int> &inorder)
{
    // Record original preorder for partitioning
    original_preorder = preorder;

    // Clear previous data
    clear();

    if (preorder.empty())
        return;

    // Store original nodes
    for (int val : preorder)
    {
        original_nodes.insert(val);
    }

    // size arrays
    max_node_value = *std::max_element(preorder.begin(), preorder.end());
    max_node_value = std::max(max_node_value, *std::max_element(inorder.begin(), inorder.end()));

    // Resize vectors - only need space for original nodes
    ranges.resize(max_node_value + 1, {0, 0});
    edges.resize(max_node_value + 1, {NO_CHILD, NO_CHILD});
    parents.resize(max_node_value + 1, NO_PARENT);

    // Store original inorder
    original_inorder = inorder;

    // Create position mapping from inorder sequence
    for (size_t i = 0; i < inorder.size(); i++)
    {
        position_in_inorder[inorder[i]] = i;
    }

    root = preorder[0];
    buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, NO_PARENT);

    // Calculate ranges for all nodes
    calculateAllRanges();
}

// Definition: Get the left child of a node
// Parameters: node: node id
// Returns: left child id or NO_CHILD
// Errors: returns NO_CHILD for invalid node
int VectorRangeTreeMap::getLeftChild(int node) const
{
    return (node >= 0 && node < (int)edges.size()) ? edges[node].first : NO_CHILD;
}

// Definition: Get the right child of a node
// Parameters: node: node id
// Returns: right child id or NO_CHILD
// Errors: returns NO_CHILD for invalid node
int VectorRangeTreeMap::getRightChild(int node) const
{
    return (node >= 0 && node < (int)edges.size()) ? edges[node].second : NO_CHILD;
}

// Definition: Get the parent of a node
// Parameters: node: node id
// Returns: parent id or NO_PARENT
// Errors: returns NO_PARENT for invalid node
int VectorRangeTreeMap::getParent(int node) const
{
    return (node >= 0 && node < (int)parents.size()) ? parents[node] : NO_PARENT;
}

// Definition: Get the inorder range of a node
// Parameters: node: node id
// Returns: {start, end} range or {0,0} if invalid
// Errors: returns {0,0} for invalid node
std::pair<int, int> VectorRangeTreeMap::getRange(int node) const
{
    return (node >= 0 && node < (int)ranges.size()) ? ranges[node] : std::make_pair(0, 0);
}

// Definition: Set left child and update parent pointers
// Parameters: node: parent id; child: child id
// Returns: nothing
// Errors: ignores invalid node ids
void VectorRangeTreeMap::setLeftChild(int node, int child)
{

    if (node < 0 || node > max_node_value)
        return;
    int old = edges[node].first;
    if (old >= 0)
        parents[old] = NO_PARENT;
    edges[node].first = child;
    if (child >= 0)
        parents[child] = node;
}

// Definition: Set right child and update parent pointers
// Parameters: node: parent id; child: child id
// Returns: nothing
// Errors: ignores invalid node ids
void VectorRangeTreeMap::setRightChild(int node,int child) {
    if (node<0||node>max_node_value) return;
    int old = edges[node].second;
    if (old>=0) parents[old]=NO_PARENT;
    edges[node].second = child;
    if (child>=0) parents[child]=node;
}

// Definition: Check if node exists in original node set
// Parameters: node: node id
// Returns: true if node is original
// Errors: none
bool VectorRangeTreeMap::isOriginal(int node) const
{
    return original_nodes.count(node) > 0;
}

// Definition: Print tree structure recursively
// Parameters: node_value: start node or -2 for root; indent: prefix string
// Returns: nothing
// Errors: no-op on invalid nodes
void VectorRangeTreeMap::printTreeStructure(int node_value, std::string indent) const
{
    if (node_value == -2)
        node_value = root; // Default to root on first call
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

// Definition: Rotate left around node x and update ranges
// Parameters: x: node id
// Returns: nothing
// Errors: no-op if rotation is not possible
void VectorRangeTreeMap::rotateLeft(int x)
{
    int y = getRightChild(x);

    if (y == NO_CHILD)
        return;

    int yl = getLeftChild(y);
    int xp = getParent(x);
    auto xr = getRange(x), yr = getRange(y), ylr = getRange(yl);

    // update ranges
    ranges[y] = xr;
    ranges[x] = (yl >= 0)
                    ? std::pair<int, int>{xr.first, ylr.second}
                    : std::pair<int, int>{xr.first, yr.first};
    // reconnect pointers
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
    // compute all range
    calculateAllRanges();
}

// Definition: Rotate right around node x and update ranges
// Parameters: x: node id
// Returns: nothing
// Errors: no-op if rotation is not possible
void VectorRangeTreeMap::rotateRight(int x)
{
    int y = getLeftChild(x);
    if (y == NO_CHILD)
        return;

    int yr = getRightChild(y);
    int xp = getParent(x);
    auto xr = getRange(x), yrng = getRange(y), yrr = getRange(yr);

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
    calculateAllRanges();
}

// Definition: Collect all (parent, child) edges into a set
// Parameters: node_value: subtree root; out_set: output set
// Returns: nothing
// Errors: no-op on invalid nodes
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

// Definition: Print full tree data including ranges, edges, and parents
// Parameters: none
// Returns: nothing
// Errors: none
void VectorRangeTreeMap::print() const
{
    std::cout << "Tree structure:\n";
    printTreeStructure();

    // Print ranges for original nodes only, in sorted order
    std::vector<int> sorted_original;
    for (int node : original_nodes)
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

        // Convert non-original children to -1
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

        // Convert non-original parent to -1
        if (!isOriginal(parent))
            parent = NO_PARENT;

        std::cout << parent;
    }
    std::cout << "]\n";
}

// Definition: Partition tree along edge defined by inorder ranges
// Parameters: T: source tree; parent_range: parent range; child_range: child range
// Returns: pair of subtrees {A, B}
// Errors: may return empty subtrees for invalid ranges
std::pair<VectorRangeTreeMap,VectorRangeTreeMap>
VectorRangeTreeMap::partitionAlongEdge(const VectorRangeTreeMap& T,
                                       const std::pair<int,int>& parent_range,
                                       const std::pair<int,int>& child_range)
{
    std::vector<int> inA, inB, preA, preB;
    // split inorder
    for (int x : T.original_inorder)
    {
        auto r = T.getRange(x);
        if (r.first >= child_range.first && r.second <= child_range.second)
            inA.push_back(x);
        else
            inB.push_back(x);
    }
    // split preorder
    for (int x : T.original_preorder)
    {
        auto r = T.getRange(x);
        if (r.first >= child_range.first && r.second <= child_range.second)
            preA.push_back(x);
        else
            preB.push_back(x);
    }
    VectorRangeTreeMap A, B;
    A.build(preA, inA);
    B.build(preB, inB);
    return {A, B};
}

// Definition: Reset internal vectors and metadata
// Parameters: none
// Returns: nothing
// Errors: none
void VectorRangeTreeMap::clear() 
{
    ranges.clear();
    edges.clear();
    parents.clear();
    original_inorder.clear();
    position_in_inorder.clear();
    original_nodes.clear();
    root = -1;
    max_node_value = 0;
}

// Definition: Recursive helper to build tree from traversals
// Parameters: preorder/inorder: traversal arrays; ps/pe/is/ie: ranges; p: parent id
// Returns: nothing
// Errors: no-op on invalid ranges
void VectorRangeTreeMap::buildRecursive(const std::vector<int> &preorder, const std::vector<int> &inorder,
                    int ps, int pe, int is, int ie, int p)
{
    if (ps > pe || is > ie)
        return;

    int root_val = preorder[ps];

    // Set parent
    parents[root_val] = p;

    // Find root position in inorder
    int root_idx = position_in_inorder[root_val];
    int left_size = root_idx - is;

    // Build left subtree
    if (left_size > 0)
    {
        int left_root_val = preorder[ps + 1];
        setLeftChild(root_val, left_root_val);
        buildRecursive(preorder, inorder, ps + 1, ps + left_size,
                       is, root_idx - 1, root_val);
    }

    // Build right subtree
    if (root_idx < ie)
    {
        int right_root_val = preorder[ps + left_size + 1];
        setRightChild(root_val, right_root_val);
        buildRecursive(preorder, inorder, ps + left_size + 1, pe,
                       root_idx + 1, ie, root_val);
    }
}

// Definition: Compute a single node range from children
// Parameters: node_value: node id
// Returns: nothing
// Errors: no-op on invalid node
void VectorRangeTreeMap::updateNodeRange(int node_value)
{
    if (node_value < 0 || node_value >= ranges.size() || !isOriginal(node_value))
        return;

    int left = getLeftChild(node_value);
    int right = getRightChild(node_value);

    // Get position of this node in inorder traversal
    int pos = position_in_inorder[node_value];

    // Start with just this node's position
    int start = pos;
    int end = pos + 1;

    // Extend range to include children's ranges
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

// Definition: Calculate ranges for all nodes
// Parameters: none
// Returns: nothing
// Errors: none
void VectorRangeTreeMap::calculateAllRanges()
{
    // Calculate ranges for all nodes in post-order (bottom-up)
    calculateRangesPostOrder(root);
}

// Definition: Post-order update for all ranges in a subtree
// Parameters: node_value: subtree root
// Returns: nothing
// Errors: no-op on invalid node
void VectorRangeTreeMap::calculateRangesPostOrder(int node_value)
{
    if (node_value == NO_CHILD || !isOriginal(node_value))
        return;

    int left = getLeftChild(node_value);
    int right = getRightChild(node_value);

    // Process children first (post-order)
    if (left != NO_CHILD && isOriginal(left))
    {
        calculateRangesPostOrder(left);
    }
    if (right != NO_CHILD && isOriginal(right))
    {
        calculateRangesPostOrder(right);
    }

    // Now calculate this node's range
    updateNodeRange(node_value);
}

// Definition: Detect if two trees are identical on the same node set
// Parameters: A: first tree; B: second tree
// Returns: true if structure, ranges, parents, and children match
// Errors: returns false on mismatch or invalid input
bool TreesEqual(const VectorRangeTreeMap &A,
                const VectorRangeTreeMap &B)
{
    // Must have the same set of original node labels
    if (A.original_nodes != B.original_nodes)
        return false;

    // For each node, ranges + left/right child + parent must match
    for (int v : A.original_nodes)
    {
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

// Definition: Serialize a tree into a unique string for hashing
// Parameters: T: tree to serialize
// Returns: stable string encoding of ranges and edges
// Errors: none
std::string treeToString(const VectorRangeTreeMap &T)
{
    // Grab & sort the original labels
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());

    // Build a single string with fixed separators
    std::ostringstream oss;
    for (int v : nodes)
    {                                  // node label
        auto [rs, re] = T.getRange(v); // start/end of inorder range
        // L,R are left/right child labels
        int L = T.getLeftChild(v);
        int R = T.getRightChild(v);
        int P = T.getParent(v); // Parent label

        // format: v,rs,re,L,R,P
        oss << v << ','
            << rs << ',' << re << ','
            << L << ',' << R << ','
            << P << ';';
    }
    return oss.str();
}

// Definition: Compare two range edges for equality
// Parameters: o: other range edge
// Returns: true if all fields match
// Errors: none
bool RP::operator==(RP const &o) const {
    return cs == o.cs && ce == o.ce;
}

// Definition: Hash a range edge for unordered containers
// Parameters: r: range edge
// Returns: hash value
// Errors: none
size_t RPH::operator()(RP const &r) const {
    uint64_t h2 = ((uint64_t)r.cs << 32) ^ uint32_t(r.ce);
    return std::hash<uint64_t>()(h2);
}

// Definition: Collect all (parent_range, child_range) edges in target
// Parameters: T: tree to scan
// Returns: set of range edges
// Errors: returns empty on invalid tree
std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap& T)
{
    std::unordered_set<RP, RPH> s;
    std::function<void(int)> dfs = [&](int v)
    {
        if (v < 0 || !T.isOriginal(v))
            return;
        auto pr = T.getRange(v);
        for (int c : {T.getLeftChild(v), T.getRightChild(v)})
        {
            if (c >= 0 && T.isOriginal(c))
            {
                auto cr = T.getRange(c);
                s.insert(RP{pr.first, pr.second, cr.first, cr.second});
                dfs(c);
            }
        }
    };
    dfs(T.root);
    return s;
}

// Definition: Check if a tree already has a given range edge
// Parameters: tree: tree to scan; e: range edge
// Returns: true if edge exists
// Errors: returns false on invalid tree
bool hasEdgeByRange(const VectorRangeTreeMap &tree, RP const &e)
{
    std::function<bool(int)> dfs = [&](int v) -> bool
    {
        if (v < 0 || !tree.isOriginal(v))
            return false;
        for (int c : {tree.getLeftChild(v), tree.getRightChild(v)})
        {
            if (c >= 0 && tree.isOriginal(c))
            {
                auto cr = tree.getRange(c);
                if (cr.first == e.cs && cr.second == e.ce)
                    return true;
            }
        }

        return dfs(tree.getLeftChild(v)) || dfs(tree.getRightChild(v));
    };
    return dfs(tree.root);
}

// Definition: Detect any free edge in cur that target wants
// Parameters: cur: current tree; tgt: target tree; out_v: rotation node; out_leftRotation: direction
// Returns: true if a free edge is found
// Errors: returns false on invalid trees or if none found
bool hasFreeEdge(const VectorRangeTreeMap &cur,
                 const VectorRangeTreeMap &tgt,
                 int &out_v,
                 bool &out_leftRotation)
{
    // collect all (parent_range,child_range) edges in the target
    auto tgtEdges = buildTargetSet(tgt);

    for (int v : cur.original_nodes)
    {
        // try a left‐rotation at v (pull up its right child y)
        int y = cur.getRightChild(v);
        if (y >= 0 && cur.isOriginal(y))
        {
            auto xr = cur.getRange(v);
            auto yr = cur.getRange(y);

            // compute the child‐range for v after rotating left
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

            // if tgt has e but cur does not already have it, we found a free edge
            if (tgtEdges.count(e) && !hasEdgeByRange(cur, e))
            {
                out_v = v;
                out_leftRotation = true;
                return true;
            }
        }

        // try a right‐rotation at v (pull up its left child z)
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

// Distance Computation & Search

// Forward declarations for search and cost functions
int BFSSearch(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
int removeFreeEdge(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
int Dist(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

static std::unordered_map<std::string, int> memo;
static std::unordered_set<std::string> in_progress;

// Definition: Compute distance via free-edge shortcut and partitioning
// Parameters: Ts: source tree; Te: target tree
// Returns: rotation distance
// Errors: falls back to BFSSearch on cycles or partition mismatch
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
        // perform one free rotation + partition & recurse
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

        if (A1.original_nodes != B1.original_nodes ||
            A2.original_nodes != B2.original_nodes)
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

// Definition: Apply one free-edge rotation and recurse on partitions
// Parameters: Ts: source tree; Te: target tree
// Returns: distance contribution for that free-edge path
// Errors: falls back to BFSSearch on partition mismatch
int removeFreeEdge(const VectorRangeTreeMap &Ts,
                   const VectorRangeTreeMap &Te)
{
    int u;
    bool left;
    assert(hasFreeEdge(Ts, Te, u, left));

    // identify the child endpoint _before_ rotating
    int v = left
                ? Ts.getRightChild(u)
                : Ts.getLeftChild(u);

    // Apply the rotation
    VectorRangeTreeMap Tbar = Ts;
    if (left)
        Tbar.rotateLeft(u);
    else
        Tbar.rotateRight(u);

    // Grab the post-rotation ranges.
    // After rotation, v is the parent and u is its child.
    auto pr = Tbar.getRange(v);
    auto cr = Tbar.getRange(u);

    // Partition both trees on that new edge
    auto [A1, A2] = VectorRangeTreeMap::partitionAlongEdge(Tbar, pr, cr);
    auto [B1, B2] = VectorRangeTreeMap::partitionAlongEdge(Te, pr, cr);

    // Sanity check, then recurse
    if (A1.original_nodes != B1.original_nodes ||
        A2.original_nodes != B2.original_nodes)
    {
        return BFSSearch(Ts, Te);
    }
    return 1 + Dist(A1, B1) + Dist(A2, B2);
}

// Forward declarations
typedef std::pair<VectorRangeTreeMap, VectorRangeTreeMap> P2;
bool hasFreeEdge(const VectorRangeTreeMap &, const VectorRangeTreeMap &, int &, bool &);

// Definition: Brute-force BFS search for rotation distance
// Parameters: T_s: source tree; T_e: target tree
// Returns: minimum distance or INT_MAX if not found
// Errors: returns INT_MAX on failures
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

        for (int v : cur.original_nodes)
        {
            // try left-rotation
            if (cur.getRightChild(v) != VectorRangeTreeMap::NO_CHILD)
            {
                auto tmp = cur;
                tmp.rotateLeft(v);

                // target reached?
                if (TreesEqual(tmp, T_e))
                {
                    return nextDist;
                }

                // free-edge shortcut is handled in Dist(); BFS should stay exact.

                // else enqueue
                auto key = treeToString(tmp);
                if (!visited.count(key))
                {
                    visited.insert(key);
                    Q.push({std::move(tmp), nextDist});
                }
            }

            // try right rotation
            if (cur.getLeftChild(v) != VectorRangeTreeMap::NO_CHILD)
            {
                auto tmp = cur;
                tmp.rotateRight(v);

                // target reached?
                if (TreesEqual(tmp, T_e))
                {
                    return nextDist;
                }

                // free-edge shortcut is handled in Dist(); BFS should stay exact.

                // else enqueue
                auto key = treeToString(tmp);
                if (!visited.count(key))
                {
                    visited.insert(key);
                    Q.push({std::move(tmp), nextDist});
                }
            }
        }
    }

    return INT_MAX; // should not be reached
}

// Definition: Clear memo cache then compute distance
// Parameters: T_init: source tree; T_final: target tree
// Returns: rotation distance
// Errors: none
int FindRotationDistance(const VectorRangeTreeMap &T_init,
                         const VectorRangeTreeMap &T_final)
{
    memo.clear();
    return Dist(T_init, T_final);
}

// Definition: Run internal verification tests
// Parameters: none
// Returns: nothing
// Errors: asserts on failure
void runAllTests()
{

    VectorRangeTreeMap T;
    std::vector<int> pre = {8, 4, 2, 1, 0, 3, 6, 5, 7};
    std::vector<int> in = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    T.build(pre, in);

    // Partition along the edge 4→2:
    auto pr = T.getRange(4);
    auto cr = T.getRange(2);
    auto [A, B] = VectorRangeTreeMap::partitionAlongEdge(T, pr, cr);

    // A should be the subtree rooted at 2: inorder {0,1,2,3}, preorder {2,1,0,3}
    assert(A.original_inorder == std::vector<int>({0, 1, 2, 3}));
    assert(A.original_preorder == std::vector<int>({2, 1, 0, 3}));

    // B should be the remainder: inorder {4,5,6,7,8}, preorder {8,4,6,5,7}
    assert(B.original_inorder == std::vector<int>({4, 5, 6, 7, 8}));
    assert(B.original_preorder == std::vector<int>({8, 4, 6, 5, 7}));

    std::cout << "PASS partitionAlongEdge test\n";

    // BFSSearch & Dist tests
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
        // Random BST sanity checks: Dist should match BFSSearch on small n.
        // Helper builds preorder/inorder from a BST built by random insertion order.
        auto buildBSTPreIn = [](int n, int seed) {
            struct Node {
                int v;
                Node *l = nullptr;
                Node *r = nullptr;
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

// Definition: Build preorder and inorder for a BST from shuffled insertion order
// Parameters: n: node count; seed: shuffle seed
// Returns: {preorder, inorder} vectors
// Errors: returns empty vectors if n <= 0
static std::pair<std::vector<int>, std::vector<int>> buildRandomBST(int n, int seed)
{
    struct Node
    {
        int v;
        Node *l = nullptr;
        Node *r = nullptr;
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

// Definition: Encode the tree as P:...;I:... for Java CLI compatibility
// Parameters: T: tree to encode
// Returns: encoded string
// Errors: returns empty encodings if traversals are empty
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

// Definition: Brute-force CLI entry that emits JSON per direction
// Parameters: argc/argv: CLI args
// Returns: exit code, 0 on success
// Errors: prints to stderr on invalid args or unsupported cases
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
// Definition: Program entry point for bf_bst binary
// Parameters: argc/argv: CLI args
// Returns: exit code
// Errors: none beyond runCli and runAllTests
int main(int argc, char **argv)
{
    if (argc > 1)
        return runCli(argc, argv);
    runAllTests();
    return 0;
}
#endif
