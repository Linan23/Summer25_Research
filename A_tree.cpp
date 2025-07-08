#include "A_tree.h"
#include <iostream>
#include <algorithm>
#include <set>
#include <sstream>
#include <queue>
#include <climits>
#include <cassert>
#include <functional>

//------------------------------------------------------------------------------
// Utility for hashing (parent, child) pairs
//   Provides hashing and equality for (parent,child) pairs so we can
//   put edges into unordered_sets quickly.
//------------------------------------------------------------------------------
/*
Helper function that takes a (parent, child) pair and turns it into a single number so the hash-set can store it quickly.
It mixes the two integers into one 64-bit value and runs the standard hash on it
*/
size_t PairHash::operator()(const std::pair<int, int> &p) const
{
    // Combine two 32-bit ints into one 64-bit for hashing
    return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned long long)p.second);
}

/*
 Helper function to show how hash-set checks if two (parent, child) pairs are exactly the same.
 It returns true only if both the parent and the child values match.
 */
bool PairEq::operator()(const std::pair<int, int> &a,
                        const std::pair<int, int> &b) const
{
    // Return true only if both parent and child match
    return a.first == b.first && a.second == b.second;
}

//------------------------------------------------------------------------------
// VectorRangeTreeMap
// Two-vector BST representation + subtree-range info for partitioning
//------------------------------------------------------------------------------

// Constructor
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

// Build from preorder + inorder
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

// Accessors

// Get left child of a node
int VectorRangeTreeMap::getLeftChild(int node) const
{
    return (node >= 0 && node < (int)edges.size()) ? edges[node].first : NO_CHILD;
}

// Get right child of a node
int VectorRangeTreeMap::getRightChild(int node) const
{
    return (node >= 0 && node < (int)edges.size()) ? edges[node].second : NO_CHILD;
}

// Get parent of a node
int VectorRangeTreeMap::getParent(int node) const
{
    return (node >= 0 && node < (int)parents.size()) ? parents[node] : NO_PARENT;
}

// Get range of a node
std::pair<int, int> VectorRangeTreeMap::getRange(int node) const
{
    return (node >= 0 && node < (int)ranges.size()) ? ranges[node] : std::make_pair(0, 0);
}

// Mutators

// Set left child and update parent pointers
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

// Set right child and update parent pointers
void VectorRangeTreeMap::setRightChild(int node,int child) {
    if (node<0||node>max_node_value) return;
    int old = edges[node].second;
    if (old>=0) parents[old]=NO_PARENT;
    edges[node].second = child;
    if (child>=0) parents[child]=node;
}

// Check if node is original (exists in our node set)
bool VectorRangeTreeMap::isOriginal(int node) const
{
    return original_nodes.count(node) > 0;
}

// debug: Print tree structure recursively
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

// BST rotations (maintaining ranges)

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

// Collect all (parent,child) edges into a set
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

// Print ranges/edges/parents
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

// Partition tree along edge with (parent_range,child_range)
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

// Reset internal vectors
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

// Build tree recursively
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

// Compute a single node’s range from children
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

// Calculate ranges for all nodes - only called during initial build
void VectorRangeTreeMap::calculateAllRanges()
{
    // Calculate ranges for all nodes in post-order (bottom-up)
    calculateRangesPostOrder(root);
}

// Post-order update for all ranges
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

//  Tree Comparison & Serialization

/*
Method: Detects if two trees are the same
Returns true iff A and B represent exactly the same tree on the same original nodes
*/
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

/*Method to serialize each tree into strings to use as hashkey in unordered_set(bfs)/map(memo)
Produce a unique string representing T’s shape+ranges on its original nodes
*/
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

//  Free-Edge Detection Helpers

bool RP::operator==(RP const &o) const {
    return ps == o.ps && pe == o.pe && cs == o.cs && ce == o.ce;
}

size_t RPH::operator()(RP const &r) const {
    uint64_t h1 = ((uint64_t)r.ps << 32) ^ uint32_t(r.pe);
    uint64_t h2 = ((uint64_t)r.cs << 32) ^ uint32_t(r.ce);
    return std::hash<uint64_t>()(h1 ^ (h2 * 0x9e3779b97f4a7c15ULL));
}

// collect all (parent_range,child_range) edges of target
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

// helper to see if tree already has that range-edge
bool hasEdgeByRange(const VectorRangeTreeMap &tree, RP const &e)
{
    std::function<bool(int)> dfs = [&](int v) -> bool
    {
        if (v < 0 || !tree.isOriginal(v))
            return false;
        auto pr = tree.getRange(v);
        if (pr.first == e.ps && pr.second == e.pe)
        {
            for (int c : {tree.getLeftChild(v), tree.getRightChild(v)})
            {
                if (c >= 0 && tree.isOriginal(c))
                {
                    auto cr = tree.getRange(c);
                    if (cr.first == e.cs && cr.second == e.ce)
                        return true;
                }
            }
        }

        return dfs(tree.getLeftChild(v)) || dfs(tree.getRightChild(v));
    };
    return dfs(tree.root);
}

// // detect any “free” edge in cur that target wants
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

// Recursive cost via free-edge shortcut + partitioning

int Dist(const VectorRangeTreeMap &Ts,
         const VectorRangeTreeMap &Te)
{
    if (TreesEqual(Ts, Te))
        return 0;

    std::string key = treeToString(Ts) + "|" + treeToString(Te);
    if (auto it = memo.find(key); it != memo.end())
        return it->second;

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

        auto pr = Tbar.getRange(u);
        auto cr = Tbar.getRange(v);

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
    return result;
}

// Apply one free-edge removal and recurse

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

    // Grab the post-rotation ranges
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

// Brute-force BFS search for rotation distance

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

                // free-edge shortcut?
                int u2;
                bool l2;
                if (hasFreeEdge(tmp, T_e, u2, l2))
                {
                    int cost = nextDist + removeFreeEdge(tmp, T_e);
                    return cost;
                }

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

                // free-edge shortcut?
                int u2;
                bool l2;
                if (hasFreeEdge(tmp, T_e, u2, l2))
                {
                    int cost = nextDist + removeFreeEdge(tmp, T_e);
                    return cost;
                }

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

/// Clears the memo cache then calls Dist()
int FindRotationDistance(const VectorRangeTreeMap &T_init,
                         const VectorRangeTreeMap &T_final)
{
    memo.clear();
    return Dist(T_init, T_final);
}

// verification tests
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
    std::cout << "All tests passed";
}

int main()
{
    runAllTests();
    return 0;
}
