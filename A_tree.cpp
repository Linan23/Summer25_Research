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

    int root;                                 // root node value
    int max_node_value;                      // maximum node value to size vectors

    // Store the original inorder sequence for range calculations
    std::vector<int> original_inorder;
    std::unordered_map<int, int> position_in_inorder; // node_value -> position in inorder
    std::unordered_set<int> original_nodes;  // set of original node values

    // Constants for "no child" (dummy nodes)
    static constexpr int NO_CHILD = -1;

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

        // Store original inorder
        original_inorder = inorder;

        // Create position mapping from inorder sequence
        for (size_t i = 0; i < inorder.size(); i++) {
            position_in_inorder[inorder[i]] = i;
        }

        root = preorder[0];
        buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1);

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

    // Get range of a node
    inline std::pair<int, int> getRange(int node_value) const {
        if (node_value < 0 || node_value >= ranges.size()) return {0, 0};
        return ranges[node_value];
    }

    // Set left child
    inline void setLeftChild(int node_value, int child_value) {
        if (node_value >= 0 && node_value < edges.size()) {
            edges[node_value].first = child_value;
        }
    }

    // Set right child
    inline void setRightChild(int node_value, int child_value) {
        if (node_value >= 0 && node_value < edges.size()) {
            edges[node_value].second = child_value;
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
        std::cout << indent << "Node " << node_value
                  << " (range: " << range.first << ", " << range.second << ")\n";

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

        // Update edges: y takes x's position, x becomes right child of y
        setLeftChild(x_value, y_right);  // x's left child becomes y's right child
        setRightChild(y, x_value);       // y's right child becomes x
    }

    // Left rotation at x using professor's algorithm
    void rotateLeft(int x_value) {
        if (x_value < 0 || x_value >= edges.size() || !isOriginal(x_value)) return;

        int y = getRightChild(x_value);  // right child of x
        if (y == NO_CHILD || !isOriginal(y)) return; // can't rotate if no right child

        // Get the left child of y (will become right child of x)
        int y_left = getLeftChild(y);

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

        // Update edges: y takes x's position, x becomes left child of y
        setRightChild(x_value, y_left);  // x's right child becomes y's left child
        setLeftChild(y, x_value);        // y's left child becomes x
    }

    // Collect all parent->child edges into `out_set` (only original nodes)
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
    }

private:
    void clear() {
        ranges.clear();
        edges.clear();
        original_inorder.clear();
        position_in_inorder.clear();
        original_nodes.clear();
        root = -1;
        max_node_value = 0;
    }

    // Build tree recursively
    void buildRecursive(const std::vector<int>& preorder, const std::vector<int>& inorder,
                        int pre_start, int pre_end, int in_start, int in_end) {
        if (pre_start > pre_end || in_start > in_end) return;

        int root_val = preorder[pre_start];

        // Find root position in inorder
        int root_idx = position_in_inorder[root_val];
        int left_size = root_idx - in_start;

        // Build left subtree
        if (left_size > 0) {
            int left_root_val = preorder[pre_start + 1];
            setLeftChild(root_val, left_root_val);
            buildRecursive(preorder, inorder, pre_start + 1, pre_start + left_size,
                           in_start, root_idx - 1);
        }

        // Build right subtree
        if (root_idx < in_end) {
            int right_root_val = preorder[pre_start + left_size + 1];
            setRightChild(root_val, right_root_val);
            buildRecursive(preorder, inorder, pre_start + left_size + 1, pre_end,
                           root_idx + 1, in_end);
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

// Compute and print shared edges between two VectorRangeTreeMap trees
void getSharedEdges(const VectorRangeTreeMap& A, const VectorRangeTreeMap& B) {
    // Collect edges of A
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> edgesA;
    A.collectEdges(A.root, edgesA);

    // Collect edges of B
    std::unordered_set<std::pair<int,int>, PairHash, PairEq> edgesB;
    B.collectEdges(B.root, edgesB);

    // Print intersection
    std::cout << "Shared edges (parent -> child):\n";
    for (auto& e : edgesA) {
        if (edgesB.count(e)) {
            std::cout << e.first << " -> " << e.second << "\n";
        }
    }
}

int main() {
    VectorRangeTreeMap tree;

    // Test with professor's example: inorder [0,1,2,3,4,5,6,7,8]
    std::cout << "Professor's Example Tree:\n";
    std::vector<int> preorder = {8, 4, 2, 1, 0, 3, 6, 5, 7};
    std::vector<int> inorder = {0, 1, 2, 3, 4, 5, 6, 7, 8};

    tree.build(preorder, inorder);
    std::cout << "Initial tree:\n";
    tree.print();
    std::cout << "\n";

    // Test right rotation at node 4
    std::cout << "Right rotation at node 4:\n";
    tree.rotateRight(4);
    tree.print();
    std::cout << "\n";

    // Rebuild original tree and test left rotation at node 4
    std::cout << "Left rotation at node 4 (rebuilding original first):\n";
    tree.build(preorder, inorder);
    tree.rotateLeft(4);
    tree.print();
    std::cout << "\n";

    // Another Left original tree and test left rotation at node 4
    std::cout << "Another left rotation at node 4:\n";
    tree.build(preorder, inorder);
    tree.rotateLeft(4);
    tree.print();
    std::cout << "\n";

    return 0;
}