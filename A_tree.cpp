#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <set>

/*
Helper function that takes a (parent, child) pair and turns it into a single number so the hash-set can store it quickly.
It mixes the two integers into one 64-bit value and runs the standard hash on it
*/
struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const {
        // Combine two 32-bit ints into one 64-bit for hashing
        return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned long long)p.second);
    }
};

/*
 Helper function to show how hash-set checks if two (parent, child) pairs are exactly the same.
 It returns true only if both the parent and the child values match.
 */
struct PairEq {
    bool operator()(const std::pair<int,int>& a, const std::pair<int,int>& b) const {
        return a.first == b.first && a.second == b.second;
    }
};

struct VectorRangeTreeMap {
    // Professor's two-vector representation
    std::vector<std::pair<int, int>> ranges;  // ranges[node_value] = (start, end) range
    std::vector<std::pair<int, int>> edges;   // edges[node_value] = (left_child, right_child)
    std::vector<int> parents;                 // parents[node_value] = parent_value

    int root;                                 // root node value
    int max_node_value;                      // maximum node value to size vectors

    // Store the original inorder sequence for range calculations
    std::vector<int> original_inorder;
    std::unordered_map<int, int> position_in_inorder; // node_value -> position in inorder
    std::unordered_set<int> original_nodes;  // set of original node values

    // Constants for "no child" and "no parent" (dummy nodes)
    static constexpr int NO_CHILD = -1;
    static constexpr int NO_PARENT = -1;

    VectorRangeTreeMap() : root(-1), max_node_value(0) {}

    // Build from preorder and inorder lists
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        // Clear previous data
        clear();

        if (preorder.empty()) return;

        // Store original nodes
        for (int val : preorder) {
            original_nodes.insert(val);
        }

        // Find maximum node value to size vectors
        max_node_value = *std::max_element(preorder.begin(), preorder.end());
        max_node_value = std::max(max_node_value, *std::max_element(inorder.begin(), inorder.end()));

        // Resize vectors - only need space for original nodes
        ranges.resize(max_node_value + 1, {0, 0});
        edges.resize(max_node_value + 1, {NO_CHILD, NO_CHILD});
        parents.resize(max_node_value + 1, NO_PARENT);

        // Store original inorder
        original_inorder = inorder;

        // Create position mapping from inorder sequence
        for (size_t i = 0; i < inorder.size(); i++) {
            position_in_inorder[inorder[i]] = i;
        }

        root = preorder[0];
        buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, NO_PARENT);

        // Calculate ranges for all nodes
        calculateAllRanges();
    }

    // Get left child of a node (-1 if dummy/none)
    inline int getLeftChild(int node_value) const {
        if (node_value < 0 || node_value >= edges.size()) return NO_CHILD;
        return edges[node_value].first;
    }

    // Get right child of a node (-1 if dummy/none)
    inline int getRightChild(int node_value) const {
        if (node_value < 0 || node_value >= edges.size()) return NO_CHILD;
        return edges[node_value].second;
    }

    // Get parent of a node (-1 if root/none)
    inline int getParent(int node_value) const {
        if (node_value < 0 || node_value >= parents.size()) return NO_PARENT;
        return parents[node_value];
    }

    // Get range of a node
    inline std::pair<int, int> getRange(int node_value) const {
        if (node_value < 0 || node_value >= ranges.size()) return {0, 0};
        return ranges[node_value];
    }

    // Set left child and update parent pointers
    inline void setLeftChild(int node_value, int child_value) {
        if (node_value >= 0 && node_value < edges.size()) {
            // Remove old parent relationship
            int old_child = edges[node_value].first;
            if (old_child != NO_CHILD && old_child >= 0 && old_child < parents.size()) {
                parents[old_child] = NO_PARENT;
            }

            // Set new child
            edges[node_value].first = child_value;

            // Set new parent relationship
            if (child_value != NO_CHILD && child_value >= 0 && child_value < parents.size()) {
                parents[child_value] = node_value;
            }
        }
    }

    // Set right child and update parent pointers
    inline void setRightChild(int node_value, int child_value) {
        if (node_value >= 0 && node_value < edges.size()) {
            // Remove old parent relationship
            int old_child = edges[node_value].second;
            if (old_child != NO_CHILD && old_child >= 0 && old_child < parents.size()) {
                parents[old_child] = NO_PARENT;
            }

            // Set new child
            edges[node_value].second = child_value;

            // Set new parent relationship
            if (child_value != NO_CHILD && child_value >= 0 && child_value < parents.size()) {
                parents[child_value] = node_value;
            }
        }
    }

    // Check if node is original (exists in our node set)
    inline bool isOriginal(int node_value) const {
        return original_nodes.count(node_value) > 0;
    }

    // Print tree structure recursively (only original nodes)
    void printTreeStructure(int node_value = -2, std::string indent = "") const {
        if (node_value == -2) node_value = root; // Default to root on first call
        if (node_value == NO_CHILD || !isOriginal(node_value)) return;

        auto range = getRange(node_value);
        int parent = getParent(node_value);
        std::cout << indent << "Node " << node_value
                  << " (range: " << range.first << ", " << range.second
                  << ", parent: " << parent << ")\n";

        int left = getLeftChild(node_value);
        if (left != NO_CHILD && isOriginal(left)) {
            printTreeStructure(left, indent + "  L-");
        }

        int right = getRightChild(node_value);
        if (right != NO_CHILD && isOriginal(right)) {
            printTreeStructure(right, indent + "  R-");
        }
    }

    // Right rotation at x using professor's algorithm
    void rotateRight(int x_value) {
        if (x_value < 0 || x_value >= edges.size() || !isOriginal(x_value)) return;

        int y = getLeftChild(x_value);  // left child of x
        if (y == NO_CHILD || !isOriginal(y)) return; // can't rotate if no left child

        // Get the right child of y (will become left child of x)
        int y_right = getRightChild(y);

        // Get parent of x (y will take x's place)
        int x_parent = getParent(x_value);

        // Store original ranges for professor's formula
        auto x_range = getRange(x_value);
        auto y_range = getRange(y);
        auto y_right_range = (y_right != NO_CHILD && isOriginal(y_right)) ?
                             getRange(y_right) : std::make_pair(0, 0);

        // Professor's range update formula:
        // 1. y takes x's original range
        ranges[y] = x_range;

        // 2. x gets: x_original_range - y_original_range + y_right_range
        // This is: (x_range - y_range) + y_right_range
        std::pair<int, int> new_x_range;

        if (y_right != NO_CHILD && isOriginal(y_right)) {
            // x gets what's left after removing y's range, plus y_right's range
            new_x_range = {y_right_range.first, x_range.second};
        } else {
            // No right child of y, so x gets the remainder after y's range
            new_x_range = {y_range.second, x_range.second};
        }
        ranges[x_value] = new_x_range;

        // Update parent relationships first, before changing edges
        // y will take x's place, so y gets x's parent
        parents[y] = x_parent;

        // x becomes y's right child
        parents[x_value] = y;

        // y_right becomes x's left child (if it exists)
        if (y_right != NO_CHILD && isOriginal(y_right)) {
            parents[y_right] = x_value;
        }

        // Update edges:
        // 1. x's left child becomes y's right child
        edges[x_value].first = y_right;

        // 2. y's right child becomes x
        edges[y].second = x_value;

        // 3. Update x's parent to point to y
        if (x_parent != NO_PARENT) {
            if (edges[x_parent].first == x_value) {
                edges[x_parent].first = y;
            } else {
                edges[x_parent].second = y;
            }
        } else {
            // x was the root, now y is the root
            root = y;
        }
    }

    // Left rotation at x using professor's algorithm
    void rotateLeft(int x_value) {
        if (x_value < 0 || x_value >= edges.size() || !isOriginal(x_value)) return;

        int y = getRightChild(x_value);  // right child of x
        if (y == NO_CHILD || !isOriginal(y)) return; // can't rotate if no right child

        // Get the left child of y (will become right child of x)
        int y_left = getLeftChild(y);

        // Get parent of x (y will take x's place)
        int x_parent = getParent(x_value);

        // Store original ranges for professor's formula
        auto x_range = getRange(x_value);
        auto y_range = getRange(y);
        auto y_left_range = (y_left != NO_CHILD && isOriginal(y_left)) ?
                            getRange(y_left) : std::make_pair(0, 0);

        // Professor's range update formula:
        // 1. y takes x's original range
        ranges[y] = x_range;

        // 2. x gets: x_original_range - y_original_range + y_left_range
        std::pair<int, int> new_x_range;

        if (y_left != NO_CHILD && isOriginal(y_left)) {
            // x gets from start to where y_left ends
            new_x_range = {x_range.first, y_left_range.second};
        } else {
            // No left child of y, so x gets up to where y starts
            new_x_range = {x_range.first, y_range.first};
        }
        ranges[x_value] = new_x_range;

        // Update parent relationships first, before changing edges
        // y will take x's place, so y gets x's parent
        parents[y] = x_parent;

        // x becomes y's left child
        parents[x_value] = y;

        // y_left becomes x's right child (if it exists)
        if (y_left != NO_CHILD && isOriginal(y_left)) {
            parents[y_left] = x_value;
        }

        // Update edges:
        // 1. x's right child becomes y's left child
        edges[x_value].second = y_left;

        // 2. y's left child becomes x
        edges[y].first = x_value;

        // 3. Update x's parent to point to y
        if (x_parent != NO_PARENT) {
            if (edges[x_parent].first == x_value) {
                edges[x_parent].first = y;
            } else {
                edges[x_parent].second = y;
            }
        } else {
            // x was the root, now y is the root
            root = y;
        }
    }

    // Collect all parent->child edges into out_set (only original nodes)
    void collectEdges(int node_value,
                      std::unordered_set<std::pair<int,int>, PairHash, PairEq>& out_set) const
    {
        if (node_value == NO_CHILD || !isOriginal(node_value)) return;

        int left = getLeftChild(node_value);
        if (left != NO_CHILD && isOriginal(left)) {
            out_set.insert({node_value, left});
            collectEdges(left, out_set);
        }

        int right = getRightChild(node_value);
        if (right != NO_CHILD && isOriginal(right)) {
            out_set.insert({node_value, right});
            collectEdges(right, out_set);
        }
    }

    // Print ranges and edges in professor's format (only original nodes)
    void print() const {
        std::cout << "Tree structure:\n";
        printTreeStructure();

        // Print ranges for original nodes only, in sorted order
        std::vector<int> sorted_original;
        for (int node : original_nodes) {
            sorted_original.push_back(node);
        }
        std::sort(sorted_original.begin(), sorted_original.end());

        std::cout << "\nranges: [";
        for (size_t i = 0; i < sorted_original.size(); i++) {
            if (i > 0) std::cout << ", ";
            auto range = getRange(sorted_original[i]);
            std::cout << "(" << range.first << "," << range.second << ")";
        }
        std::cout << "]\n";

        std::cout << "edges: [";
        for (size_t i = 0; i < sorted_original.size(); i++) {
            if (i > 0) std::cout << ", ";
            int node = sorted_original[i];
            int left = getLeftChild(node);
            int right = getRightChild(node);

            // Convert non-original children to -1
            if (!isOriginal(left)) left = NO_CHILD;
            if (!isOriginal(right)) right = NO_CHILD;

            std::cout << "(" << left << "," << right << ")";
        }
        std::cout << "]\n";

        std::cout << "parents: [";
        for (size_t i = 0; i < sorted_original.size(); i++) {
            if (i > 0) std::cout << ", ";
            int node = sorted_original[i];
            int parent = getParent(node);

            // Convert non-original parent to -1
            if (!isOriginal(parent)) parent = NO_PARENT;

            std::cout << parent;
        }
        std::cout << "]\n";
    }

private:
    void clear() {
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
    void buildRecursive(const std::vector<int>& preorder, const std::vector<int>& inorder,
                        int pre_start, int pre_end, int in_start, int in_end, int parent_val) {
        if (pre_start > pre_end || in_start > in_end) return;

        int root_val = preorder[pre_start];

        // Set parent
        parents[root_val] = parent_val;

        // Find root position in inorder
        int root_idx = position_in_inorder[root_val];
        int left_size = root_idx - in_start;

        // Build left subtree
        if (left_size > 0) {
            int left_root_val = preorder[pre_start + 1];
            setLeftChild(root_val, left_root_val);
            buildRecursive(preorder, inorder, pre_start + 1, pre_start + left_size,
                           in_start, root_idx - 1, root_val);
        }

        // Build right subtree
        if (root_idx < in_end) {
            int right_root_val = preorder[pre_start + left_size + 1];
            setRightChild(root_val, right_root_val);
            buildRecursive(preorder, inorder, pre_start + left_size + 1, pre_end,
                           root_idx + 1, in_end, root_val);
        }
    }

    // Calculate range for a single node based on its subtree span in inorder
    void updateNodeRange(int node_value) {
        if (node_value < 0 || node_value >= ranges.size() || !isOriginal(node_value)) return;

        int left = getLeftChild(node_value);
        int right = getRightChild(node_value);

        // Get position of this node in inorder traversal
        int pos = position_in_inorder[node_value];

        // Start with just this node's position
        int start = pos;
        int end = pos + 1;

        // Extend range to include children's ranges
        if (left != NO_CHILD && isOriginal(left)) {
            auto left_range = getRange(left);
            start = left_range.first;
        }

        if (right != NO_CHILD && isOriginal(right)) {
            auto right_range = getRange(right);
            end = right_range.second;
        }

        ranges[node_value] = {start, end};
    }

    // Calculate ranges for all nodes - only called during initial build
    void calculateAllRanges() {
        // Calculate ranges for all nodes in post-order (bottom-up)
        calculateRangesPostOrder(root);
    }

    // Post-order traversal to calculate ranges bottom-up
    void calculateRangesPostOrder(int node_value) {
        if (node_value == NO_CHILD || !isOriginal(node_value)) return;

        int left = getLeftChild(node_value);
        int right = getRightChild(node_value);

        // Process children first (post-order)
        if (left != NO_CHILD && isOriginal(left)) {
            calculateRangesPostOrder(left);
        }
        if (right != NO_CHILD && isOriginal(right)) {
            calculateRangesPostOrder(right);
        }

        // Now calculate this node's range
        updateNodeRange(node_value);
    }
};

/*
// Comparison function Compute and print shared edges between two VectorRangeTreeMap trees 
void getSharedEdges(const VectorRangeTreeMap& A, const VectorRangeTreeMap& B) {
    // Collect all parent→child pairs from tree A into a hash‐set 
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> edgesA;
    A.collectEdges(A.root, edgesA);  // uses collectEdges implementation

    // Same as previous
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> edgesB;
    B.collectEdges(B.root, edgesB);  // uses collectEdges implementation

    // For each (parent,child) in edgesA, check if it also appears in edgesB
    std::cout << "Shared edges (parent -> child):\n";
    for (auto& e : edgesA) {
        if (edgesB.count(e)) {
            std::cout << e.first << " -> " << e.second << "\n";
        }
    }
}
*/

// “By‐range” comparison instead of label‐comparison
void getSharedEdgesByRange(
    const VectorRangeTreeMap& A,
    const VectorRangeTreeMap& B)
{
    // 1) Build a map in Tree A from (parentRange, childRange) → (parentLabel, childLabel)
    //
    //    We pack two intervals ((ps,pe),(cs,ce)) into a single key struct.
    struct RangePair {
        int ps, pe;   // parent range [start,end)
        int cs, ce;   // child  range [start,end)
        bool operator==(const RangePair& o) const {
            return ps==o.ps && pe==o.pe
                && cs==o.cs && ce==o.ce;
        }
    };

    struct RangePairHash {
        std::size_t operator()(const RangePair& r) const {
            // Combine four 32‐bit ints into two 64‐bit values, then mix
            uint64_t h1 = ((uint64_t)r.ps << 32) ^ (uint32_t)r.pe;
            uint64_t h2 = ((uint64_t)r.cs << 32) ^ (uint32_t)r.ce;
            // A basic “xor & std::hash” mix:
            return std::hash<uint64_t>()(h1 ^ (h2 * 0x9e3779b97f4a7c15ULL));
        }
    };

    // Map from RangePair → (parentLabel, childLabel) for Tree A
    std::unordered_map<RangePair, std::pair<int,int>, RangePairHash> mapA;
    mapA.reserve( A.original_nodes.size() * 2 );

    // Recursively collect all edges p→c in A, storing their ranges:
    std::function<void(int)> collectA = [&](int node) {
        if (node < 0 || !A.isOriginal(node)) return;

        auto pr = A.getRange(node);
        int ps = pr.first, pe = pr.second;

        int L = A.getLeftChild(node);
        if (L != VectorRangeTreeMap::NO_CHILD && A.isOriginal(L)) {
            auto cr = A.getRange(L);
            RangePair key{ ps, pe, cr.first, cr.second };
            mapA[key] = { node, L };
            collectA(L);
        }

        int R = A.getRightChild(node);
        if (R != VectorRangeTreeMap::NO_CHILD && A.isOriginal(R)) {
            auto cr = A.getRange(R);
            RangePair key{ ps, pe, cr.first, cr.second };
            mapA[key] = { node, R };
            collectA(R);
        }
    };
    if (A.root != VectorRangeTreeMap::NO_CHILD) {
        collectA(A.root);
    }

    // 2) Walk Tree B. For each edge q→d, form its (q.range,d.range) and look in mapA.
    std::cout << "Shared edges (parent -> child) by RANGE:\n";

    std::function<void(int)> collectB = [&](int node) {
        if (node < 0 || !B.isOriginal(node)) return;

        auto pr = B.getRange(node);
        int ps = pr.first, pe = pr.second;

        int L = B.getLeftChild(node);
        if (L != VectorRangeTreeMap::NO_CHILD && B.isOriginal(L)) {
            auto cr = B.getRange(L);
            RangePair key{ ps, pe, cr.first, cr.second };
            auto it = mapA.find(key);
            if (it != mapA.end()) {
                auto [pa, ca] = it->second;
                std::cout << "A’s (" << pa << " -> " << ca << ")"
                          << "  matches  "
                          << "B’s (" << node << " -> " << L << ")\n";
            }
            collectB(L);
        }

        int R = B.getRightChild(node);
        if (R != VectorRangeTreeMap::NO_CHILD && B.isOriginal(R)) {
            auto cr = B.getRange(R);
            RangePair key{ ps, pe, cr.first, cr.second };
            auto it = mapA.find(key);
            if (it != mapA.end()) {
                auto [pa, ca] = it->second;
                std::cout << "A’s (" << pa << " -> " << ca << ")"
                          << "  matches  "
                          << "B’s (" << node << " -> " << R << ")\n";
            }
            collectB(R);
        }
    };
    if (B.root != VectorRangeTreeMap::NO_CHILD) {
        collectB(B.root);
    }
}


int main() {
    VectorRangeTreeMap tree;

    // Test with professor's example: inorder [0,1,2,3,4,5,6,7,8]
    std::cout << "Professor's Example Tree:\n";
    std::vector<int> preorder = {8, 4, 2, 1, 0, 3, 6, 5, 7};
    std::vector<int> inorder = {0, 1, 2, 3, 4, 5, 6, 7, 8};

    tree.build(preorder, inorder);  // buildRecursive is implemented above
    std::cout << "Initial tree:\n";
    tree.print();  // print uses printTreeStructure and updateNodeRange
    std::cout << "\n";

    // Test right rotation at node 4
    std::cout << "Right rotation at node 4:\n";
    tree.rotateRight(4);  // uses rotateRight implementation
    tree.print();  
    std::cout << "\n";

    // Rebuild original tree and test left rotation at node 4
    std::cout << "Left rotation at node 4 (rebuilding original first):\n";
    tree.build(preorder, inorder);
    tree.rotateLeft(4);   // uses rotateLeft implementation
    tree.print();  
    std::cout << "\n";

    // --- Quick sanity check between two small trees ---
    std::cout << "Test: Shared edges between two trees\n";
    VectorRangeTreeMap tree1, tree2;
    tree1.build({4,2,6,5,7}, {2,4,5,6,7});
    tree2.build({4,2,6,3,5}, {2,4,3,6,5});
    std::cout << "Tree1 structure:\n"; tree1.print(); std::cout << "\n";
    std::cout << "Tree2 structure:\n"; tree2.print(); std::cout << "\n";
    getSharedEdgesByRange(tree1, tree2);

    std::cout << "Extra Unit tests for Comparison Method\n";

    {
        std::cout << "\n Case 1: Both trees empty \n";
        VectorRangeTreeMap A, B;
        A.build({}, {});
        B.build({}, {});
        std::cout << "Tree A and B are both empty.\n";
        std::cout << "Shared edges (expected none):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        std::cout << "\n Case 2: One empty, one nonempty \n";
        VectorRangeTreeMap A, B;
        A.build({}, {});                  // A is empty
        B.build({1,2,3}, {2,1,3});        // B:  1
                                          //     / \
                                          //    2   3
        std::cout << "Tree A (empty), Tree B:\n";
        B.print();
        std::cout << "Shared edges (expected none):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        std::cout << "\n Case 3: Single-node trees \n";
        VectorRangeTreeMap A, B;
        A.build({42}, {42});    // only node 42
        B.build({42}, {42});    // same single node
        std::cout << "Tree A and Tree B (single node 42):\n";
        A.print();
        std::cout << "Shared edges (expected none):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        // false postive
        std::cout << "\n Case 4: Disjoint value sets \n";
        VectorRangeTreeMap A, B;
        A.build({1,2,3}, {1,2,3});  // values {1,2,3}
        B.build({4,5,6}, {4,5,6});  // values {4,5,6}
        std::cout << "Tree A uses {1,2,3}, Tree B uses {4,5,6}.\n";
        A.print(); std::cout << "\n"; 
        B.print(); std::cout << "\n"; 
        std::cout << "Shared edges (expected none):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        std::cout << "\n Case 5: Same labels, different shape \n";
        VectorRangeTreeMap A, B;
        // A: 1→2→3→4 (right‐leaning chain)
        A.build({1,2,3,4}, {1,2,3,4});
        // B:     2
        //       / \
        //      1   3
        //           \
        //            4
        B.build({2,1,3,4}, {1,2,3,4});
        std::cout << "Tree A: right-leaning chain 1 - 2 - 3 - 4\n";
        A.print(); std::cout << "\n";
        std::cout << "Tree B: different shape, same labels\n";
        B.print(); std::cout << "\n";
        std::cout << "Shared edges (expected only “3-4”):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        std::cout << "\n Case 6: Partial overlap deep in one subtree \n";
        VectorRangeTreeMap A, B;
        // A:
        //     5
        //    / \
        //   3   8
        //  /   / \
        // 1   7   9
        A.build({5,3,1,8,7,9}, {1,3,5,7,8,9});
        // B:
        //    5
        //   / \
        //  2   8
        //     /
        //    7
        B.build({5,2,8,7}, {2,5,7,8});
        std::cout << "Tree A:\n"; A.print(); std::cout << "\n";
        std::cout << "Tree B:\n"; B.print(); std::cout << "\n";
        std::cout << "Shared edges (expected none):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        std::cout << "\n Case 7: Completely identical trees \n";
        VectorRangeTreeMap A, B;
        std::vector<int> Pre = {8,4,2,1,0,3,6,5,7};
        std::vector<int> In  = {0,1,2,3,4,5,6,7,8};
        A.build(Pre, In);
        B.build(Pre, In);
        std::cout << "Tree A and B are identical:\n";
        A.print(); std::cout << "\n";
        std::cout << "Shared edges (expected all 8 parent→child pairs):\n";
        getSharedEdgesByRange(A, B);
    }

    {
        std::cout << "\n Case 8: Larger random-like trees, small shared subtree \n";
        VectorRangeTreeMap A, B;
        // A:         10
        //           /  \
        //          5    15
        //         / \   /
        //        2   7 12
        A.build({10,5,2,4,7,15,12}, {2,4,5,7,10,12,15});
        // B:         7
        //           /  \
        //          5    10
        //         /    /  \
        //        2    -   15
        //         \        \
        //          4        12
        B.build({7,5,2,4,10,15,12}, {2,4,5,7,10,12,15});
        std::cout << "Tree A:\n"; A.print(); std::cout << "\n";
        std::cout << "Tree B:\n"; B.print(); std::cout << "\n";
        std::cout << "Shared edges (expected “2-4” and “15-12”):\n";
        getSharedEdgesByRange(A, B);
    }

    return 0;
}
