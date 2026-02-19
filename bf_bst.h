#ifndef BF_BST_H
#define BF_BST_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>

// Definition: Hash a (parent, child) pair for unordered containers
// Parameters: p: pair to hash
// Returns: hash value
// Errors: none
struct PairHash
{
    size_t operator()(const std::pair<int, int> &p) const;
};

// Definition: Compare two (parent, child) pairs for equality
// Parameters: a: first pair; b: second pair
// Returns: true if both pairs match
// Errors: none
struct PairEq
{
    bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const;
};

// Definition: Binary-tree storage with inorder range metadata
// Parameters: none
// Returns: struct instance
// Errors: none
struct VectorRangeTreeMap
{
    std::vector<std::pair<int, int>> ranges;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> parents;
    int root;
    int max_node_value;
    std::vector<int> original_inorder, original_preorder;
    std::unordered_map<int, int> position_in_inorder;
    std::unordered_set<int> original_nodes;

    static constexpr int NO_CHILD = -1;
    static constexpr int NO_PARENT = -1;

    // Definition: Initialize an empty tree with default values
    // Parameters: none
    // Returns: constructed tree
    // Errors: none
    VectorRangeTreeMap();

    // Definition: Build tree and ranges from preorder and inorder traversals
    // Parameters: preorder: preorder node list; inorder: inorder node list
    // Returns: nothing
    // Errors: leaves tree empty on invalid input
    void build(const std::vector<int> &preorder, const std::vector<int> &inorder);

public:
    // Definition: Get the left child of a node
    // Parameters: node_value: node id
    // Returns: left child id or NO_CHILD
    // Errors: returns NO_CHILD for invalid node
    int getLeftChild(int node_value) const;

    // Definition: Get the right child of a node
    // Parameters: node_value: node id
    // Returns: right child id or NO_CHILD
    // Errors: returns NO_CHILD for invalid node
    int getRightChild(int node_value) const;

    // Definition: Get the parent of a node
    // Parameters: node_value: node id
    // Returns: parent id or NO_PARENT
    // Errors: returns NO_PARENT for invalid node
    int getParent(int node_value) const;

    // Definition: Get the inorder range for a node
    // Parameters: node_value: node id
    // Returns: {start, end} range or {0,0} if invalid
    // Errors: returns {0,0} for invalid node
    std::pair<int, int> getRange(int node_value) const;

    // Definition: Check if a node is part of the original node set
    // Parameters: node_value: node id
    // Returns: true if node is original
    // Errors: none
    bool isOriginal(int node_value) const;

    // Definition: Set left child and update parent pointers
    // Parameters: node_value: parent node id; child_value: child node id
    // Returns: nothing
    // Errors: ignores invalid node ids
    void setLeftChild(int node_value, int child_value);

    // Definition: Set right child and update parent pointers
    // Parameters: node_value: parent node id; child_value: child node id
    // Returns: nothing
    // Errors: ignores invalid node ids
    void setRightChild(int node_value, int child_value);

    // Definition: Rotate a node left and update ranges
    // Parameters: x: node id to rotate
    // Returns: nothing
    // Errors: no-op if rotation is not possible
    void rotateLeft(int x);

    // Definition: Rotate a node right and update ranges
    // Parameters: x: node id to rotate
    // Returns: nothing
    // Errors: no-op if rotation is not possible
    void rotateRight(int x);

    // Definition: Partition a tree along an edge defined by inorder ranges
    // Parameters: T: source tree; parent_range: parent range; child_range: child range
    // Returns: pair of subtrees {A, B}
    // Errors: may return empty subtrees for invalid ranges
    static std::pair<VectorRangeTreeMap, VectorRangeTreeMap>
    partitionAlongEdge(const VectorRangeTreeMap &T,
                       const std::pair<int, int> &parent_range,
                       const std::pair<int, int> &child_range);

    // Definition: Collect all edges into a set from a subtree
    // Parameters: node: subtree root; out: output set
    // Returns: nothing
    // Errors: no-op on invalid nodes
    void collectEdges(int node,
                      std::unordered_set<std::pair<int,int>,PairHash,PairEq>& out) const;

    // Definition: Print the tree structure with ranges and parents
    // Parameters: node_value: starting node or -2 for root; indent: prefix string
    // Returns: nothing
    // Errors: no-op on invalid nodes
    void printTreeStructure(int node_value = -2, std::string indent = "") const;

    // Definition: Print full tree data (structure, ranges, edges, parents)
    // Parameters: none
    // Returns: nothing
    // Errors: none
    void print() const;

    // For external use
    friend bool TreesEqual(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B);
    friend std::string treeToString(const VectorRangeTreeMap &T);
    friend int BFSSearch(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
    friend int removeFreeEdge(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
    friend int Dist(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);
    friend int FindRotationDistance(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);

private:
    // Definition: Recursive builder helper for build()
    // Parameters: preorder/inorder: traversal arrays; ps/pe/is/ie: ranges; parent: parent id
    // Returns: nothing
    // Errors: no-op on invalid ranges
    void buildRecursive(const std::vector<int> &preorder, const std::vector<int> &inorder,
                        int ps, int pe, int is, int ie, int parent);

    // Definition: Reset internal storage
    // Parameters: none
    // Returns: nothing
    // Errors: none
    void clear();

    // Definition: Recompute ranges for all nodes
    // Parameters: none
    // Returns: nothing
    // Errors: none
    void calculateAllRanges();

    // Definition: Recompute range for a single node
    // Parameters: node_value: node id
    // Returns: nothing
    // Errors: no-op on invalid node
    void updateNodeRange(int node_value);

    // Definition: Post-order range update for a subtree
    // Parameters: node_value: subtree root id
    // Returns: nothing
    // Errors: no-op on invalid node
    void calculateRangesPostOrder(int node_value);
};

// Definition: Range edge key for diagonal detection (child range only)
// Parameters: none
// Returns: struct instance
// Errors: none
struct RP { int ps, pe, cs, ce; bool operator==(RP const&) const; };

// Definition: Hash for RP range edges
// Parameters: none
// Returns: struct instance
// Errors: none
struct RPH { size_t operator()(RP const&) const; };

// Definition: Collect all target range edges from a tree
// Parameters: T: tree to scan
// Returns: set of range edges
// Errors: returns empty on invalid tree
std::unordered_set<RP,RPH> buildTargetSet(const VectorRangeTreeMap&);

// Definition: Check if a tree contains a given range edge
// Parameters: tree: tree to scan; e: range edge
// Returns: true if edge exists
// Errors: returns false on invalid input
bool hasEdgeByRange(const VectorRangeTreeMap &tree, const RP &e);

// Definition: Detect a free edge that would insert a target edge
// Parameters: cur: current tree; tgt: target tree; out_v: rotation node; out_leftRotation: direction
// Returns: true if a free edge is found
// Errors: returns false on invalid trees or if no free edge exists
bool hasFreeEdge(const VectorRangeTreeMap &cur, const VectorRangeTreeMap &tgt,
                 int &out_v, bool &out_leftRotation);

// Definition: Brute-force BFS rotation distance between two trees
// Parameters: Ts: source tree; Te: target tree
// Returns: minimum rotation distance or INT_MAX if not found
// Errors: returns INT_MAX on failures
int BFSSearch(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

// Definition: Apply one free-edge rotation and recurse via partitioning
// Parameters: Ts: source tree; Te: target tree
// Returns: distance contribution for that free edge path
// Errors: falls back to BFSSearch on partition mismatch
int removeFreeEdge(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

// Definition: Compute rotation distance using free-edge shortcut + BFS fallback
// Parameters: Ts: source tree; Te: target tree
// Returns: rotation distance
// Errors: falls back to BFSSearch on cycles or mismatches
int Dist(const VectorRangeTreeMap &Ts, const VectorRangeTreeMap &Te);

// Definition: Clear memoization and compute rotation distance
// Parameters: T_init: source tree; T_final: target tree
// Returns: rotation distance
// Errors: none
int FindRotationDistance(const VectorRangeTreeMap &T_init, const VectorRangeTreeMap &T_final);

// Definition: Run internal verification tests
// Parameters: none
// Returns: nothing
// Errors: asserts on failure
void runAllTests();

#endif // BF_BST_H
