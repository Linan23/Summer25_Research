#ifndef BF_BST_H
#define BF_BST_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>
#include <cstdint>

// Purpose: Hash a (parent, child) pair for unordered containers
// Inputs: p: pair to hash
// Output: hash value
struct PairHash
{
    size_t operator()(const std::pair<int, int> &p) const;
};

// Purpose: Compare two (parent, child) pairs for equality
// Inputs: a: first pair; b: second pair
// Output: true if both pairs match
struct PairEq
{
    bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const;
};

// Purpose: Binary-tree storage with inorder range metadata
// Inputs: none
// Output: struct instance
struct VectorRangeTreeMap
{
    std::vector<std::pair<int, int>> ranges;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> parents;
    int root;
    int max_node_value;
    std::vector<int> original_inorder, original_preorder;
    std::unordered_map<int, int> position_in_inorder;
    std::vector<int> position_in_inorder_fast;
    std::unordered_set<int> original_nodes;
    std::vector<unsigned char> original_mask;
    mutable bool fingerprint_valid = false;
    mutable std::uint64_t fingerprint_h1 = 0;
    mutable std::uint64_t fingerprint_h2 = 0;

    static constexpr int NO_CHILD = -1;
    static constexpr int NO_PARENT = -1;

    // Purpose: Initialize an empty tree with default values
    // Inputs: none
    // Output: constructed tree
    VectorRangeTreeMap();

    // Purpose: Build tree and ranges from preorder and inorder traversals
    // Inputs: preorder: preorder node list; inorder: inorder node list
    // Output: nothing
    // Failure: leaves tree empty on invalid input
    void build(const std::vector<int> &preorder, const std::vector<int> &inorder);

public:
    // Purpose: Get the left child of a node
    // Inputs: node_value: node id
    // Output: left child id or NO_CHILD
    // Failure: returns NO_CHILD for invalid node
    int getLeftChild(int node_value) const;

    // Purpose: Get the right child of a node
    // Inputs: node_value: node id
    // Output: right child id or NO_CHILD
    // Failure: returns NO_CHILD for invalid node
    int getRightChild(int node_value) const;

    // Purpose: Get the parent of a node
    // Inputs: node_value: node id
    // Output: parent id or NO_PARENT
    // Failure: returns NO_PARENT for invalid node
    int getParent(int node_value) const;

    // Purpose: Get the inorder range for a node
    // Inputs: node_value: node id
    // Output: {start, end} range or {0,0} if invalid
    // Failure: returns {0,0} for invalid node
    std::pair<int, int> getRange(int node_value) const;

    // Purpose: Check if a node is part of the original node set
    // Inputs: node_value: node id
    // Output: true if node is original
    bool isOriginal(int node_value) const;

    // Purpose: Set left child and update parent pointers
    // Inputs: node_value: parent node id; child_value: child node id
    // Output: nothing
    // Failure: ignores invalid node ids
    void setLeftChild(int node_value, int child_value);

    // Purpose: Set right child and update parent pointers
    // Inputs: node_value: parent node id; child_value: child node id
    // Output: nothing
    // Failure: ignores invalid node ids
    void setRightChild(int node_value, int child_value);

    // Purpose: Rotate a node left and update ranges
    // Inputs: x: node id to rotate
    // Output: nothing
    // Failure: no-op if rotation is not possible
    void rotateLeft(int x);

    // Purpose: Rotate a node right and update ranges
    // Inputs: x: node id to rotate
    // Output: nothing
    // Failure: no-op if rotation is not possible
    void rotateRight(int x);

    // Purpose: Partition a tree along an edge defined by inorder ranges
    // Inputs: T: source tree; parent_range: parent range; child_range: child range
    // Output: pair of subtrees {A, B}
    // Failure: may return empty subtrees for invalid ranges
    static std::pair<VectorRangeTreeMap, VectorRangeTreeMap>
    partitionAlongEdge(const VectorRangeTreeMap &T,
                       const std::pair<int, int> &parent_range,
                       const std::pair<int, int> &child_range);

    // Purpose: Partition a tree when child_range is known to be a subtree range
    // Inputs: T: source tree; child_range: subtree inorder range
    // Output: pair of subtrees {A, B}
    // Failure: may return empty subtrees for invalid ranges
    static std::pair<VectorRangeTreeMap, VectorRangeTreeMap>
    partitionAlongSubtreeRange(const VectorRangeTreeMap &T,
                               const std::pair<int, int> &child_range);

    // Purpose: Collect all edges into a set from a subtree
    // Inputs: node: subtree root; out: output set
    // Output: nothing
    // Failure: no-op on invalid nodes
    void collectEdges(int node,
                      std::unordered_set<std::pair<int,int>,PairHash,PairEq>& out) const;

    // Purpose: Print the tree structure with ranges and parents
    // Inputs: node_value: starting node or -2 for root; indent: prefix string
    // Output: nothing
    // Failure: no-op on invalid nodes
    void printTreeStructure(int node_value = -2, std::string indent = "") const;

    // Purpose: Print full tree data (structure, ranges, edges, parents)
    // Inputs: none
    // Output: nothing
    void print() const;

    friend bool TreesEqual(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B);
    friend std::string treeToString(const VectorRangeTreeMap &T);
    friend int BFSSearch(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
    friend int removeFreeEdge(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
    friend int Dist(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
    friend int FindRotationDistance(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);

private:
    // Purpose: Recursive builder helper for build()
    // Inputs: preorder/inorder: traversal arrays; ps/pe/is/ie: ranges; parent: parent id
    // Output: nothing
    // Failure: no-op on invalid ranges
    void buildRecursive(const std::vector<int> &preorder, const std::vector<int> &inorder,
                        int ps, int pe, int is, int ie, int parent);

    // Purpose: Reset internal storage
    // Inputs: none
    // Output: nothing
    void clear();

    // Purpose: Shared implementation for full/light tree construction
    // Inputs: preorder/inorder traversals; populate_nodes/preorder: metadata fill toggles
    // Output: nothing
    // Failure: leaves tree empty on invalid input
    void buildImpl(const std::vector<int> &preorder, const std::vector<int> &inorder,
                   bool populate_nodes, bool store_preorder);

    // Purpose: Build tree metadata without populating original_nodes
    // Inputs: preorder/inorder: traversal arrays
    // Output: nothing
    // Failure: leaves tree empty on invalid input
    void buildLight(const std::vector<int> &preorder, const std::vector<int> &inorder);

    // Purpose: Recompute ranges for all nodes
    // Inputs: none
    // Output: nothing
    void calculateAllRanges();

    // Purpose: Recompute range for a single node
    // Inputs: node_value: node id
    // Output: nothing
    // Failure: no-op on invalid node
    void updateNodeRange(int node_value);

    // Purpose: Post-order range update for a subtree
    // Inputs: node_value: subtree root id
    // Output: nothing
    // Failure: no-op on invalid node
    void calculateRangesPostOrder(int node_value);
};

// Purpose: Range edge key for diagonal detection (child range only)
// Inputs: none
// Output: struct instance
struct RP { int ps, pe, cs, ce; bool operator==(RP const&) const; };

// Purpose: Hash for RP range edges
// Inputs: none
// Output: struct instance
struct RPH { size_t operator()(RP const&) const; };

// Purpose: Collect all target range edges from a tree
// Inputs: T: tree to scan
// Output: set of range edges
// Failure: returns empty on invalid tree
std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap&);

// Purpose: Check if a tree contains a given range edge
// Inputs: tree: tree to scan; e: range edge
// Output: true if edge exists
// Failure: returns false on invalid input
bool hasEdgeByRange(const VectorRangeTreeMap &tree, const RP &e);

// Purpose: Detect a free edge that would insert a target edge
// Inputs: cur: current tree; tgt: target tree; out_v: rotation node; out_leftRotation: direction
// Output: true if a free edge is found
// Failure: returns false on invalid trees or if no free edge exists
bool hasFreeEdge(const VectorRangeTreeMap &cur, const VectorRangeTreeMap &tgt,
                 int &out_v, bool &out_leftRotation);

// Purpose: Brute-force BFS rotation distance between two trees
// Inputs: Ts: source tree; Te: target tree
// Output: minimum rotation distance or INT_MAX if not found
// Failure: returns INT_MAX on failures
int BFSSearch(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

// Purpose: Apply one free-edge rotation and recurse via partitioning
// Inputs: Ts: source tree; Te: target tree
// Output: distance contribution for that free edge path
// Failure: falls back to BFSSearch on partition mismatch
int removeFreeEdge(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

// Purpose: Compute rotation distance using free-edge shortcut + BFS fallback
// Inputs: Ts: source tree; Te: target tree
// Output: rotation distance
// Failure: falls back to BFSSearch on cycles or mismatches
int Dist(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

// Purpose: Clear memoization and compute rotation distance
// Inputs: T_init: source tree; T_final: target tree
// Output: rotation distance
int FindRotationDistance(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);

// Purpose: Run internal verification tests
// Inputs: none
// Output: nothing
// Failure: asserts on failure
void runAllTests();

#endif
