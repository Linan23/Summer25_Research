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

struct RangeTreeMap {
    std::unordered_map<int, std::pair<int, int>> ranges; // node_value -> (start, end) range

    // Parent-child relationships (using maps instead of Node pointers)
    std::unordered_map<int, int> parent;       // child -> parent
    std::unordered_map<int, int> left_child;   // parent -> left child
    std::unordered_map<int, int> right_child;  // parent -> right child
    int root;                                  // root node value
    std::vector<int> leaves;                   // leaf values in left-to-right order
    std::unordered_map<int, int> nodes;       // map value -> node (for compatibility with pointer version)

    // Store the original inorder sequence for range calculations
    std::vector<int> original_inorder;

    // Set to track which nodes are dummy leaves (negative values)
    std::unordered_set<int> dummy_leaves;

    // Position mapping for O(1) range calculations
    std::unordered_map<int, int> position_map;

    // Next dummy leaf ID (negative numbers)
    int next_dummy_id;

    RangeTreeMap() : root(-1), next_dummy_id(-1000) {}

    // Build from preorder and inorder lists
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        // Clear previous data (leaves and nodes)
        ranges.clear();
        parent.clear();
        left_child.clear();
        right_child.clear();
        dummy_leaves.clear();
        position_map.clear();
        leaves.clear();
        nodes.clear();
        root = -1;
        next_dummy_id = -1000;

        // Store the original inorder sequence
        original_inorder = inorder;

        // Create position mapping from ORIGINAL inorder sequence
        for (size_t i = 0; i < original_inorder.size(); i++) {
            position_map[original_inorder[i]] = i;
        }

        if (!preorder.empty()) {
            root = preorder[0];
            buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, -1);

            // Build nodes map for compatibility
            for (int node : original_inorder) {
                nodes[node] = node; // In map version, we just store the value itself
            }

            // Build leaves vector in left-to-right order
            buildLeavesVector(root);

            addDummyLeaves();
            calculateAllRanges();
        }
    }

    // Print each node and its range
    void printTreeStructure(int node = -2, std::string indent = "") const {
        if (node == -2) node = root; // Default to root on first call
        if (node == -1) return;

        auto range_it = ranges.find(node);
        if (range_it != ranges.end()) {
            std::cout << indent
                      << "Node " << node
                      << " (range: " << range_it->second.first
                      << ", " << range_it->second.second << ")\n";
        }

        auto left_it = left_child.find(node);
        if (left_it != left_child.end() && dummy_leaves.find(left_it->second) == dummy_leaves.end()) {
            printTreeStructure(left_it->second, indent + "  L-");
        }

        auto right_it = right_child.find(node);
        if (right_it != right_child.end() && dummy_leaves.find(right_it->second) == dummy_leaves.end()) {
            printTreeStructure(right_it->second, indent + "  R-");
        }
    }

    // Print just the leaf values in order
    void printLeaves() const {
        std::cout << "Leaves (left to right): ";
        for(int leaf: leaves)
            std::cout << leaf << " ";
        std::cout << "\n";
    }

    // Rotate x and its right child up (left-rotate)
    void rotateLeft(int x) {
        if (x == -1) return;

        auto right_it = right_child.find(x);
        if (right_it == right_child.end()) return; // no right child

        int y = right_it->second;
        if (dummy_leaves.find(y) != dummy_leaves.end()) return; // can't rotate dummy

        // Perform rotation logic
        auto y_left_it = left_child.find(y);
        int y_left = (y_left_it != left_child.end()) ? y_left_it->second : -1;

        // Update x's right child
        if (y_left != -1) {
            right_child[x] = y_left;
            parent[y_left] = x;
        } else {
            right_child.erase(x);
        }

        // Update y's parent
        auto x_parent_it = parent.find(x);
        if (x_parent_it != parent.end()) {
            int x_parent = x_parent_it->second;
            parent[y] = x_parent;

            auto x_parent_left_it = left_child.find(x_parent);
            if (x_parent_left_it != left_child.end() && x_parent_left_it->second == x) {
                left_child[x_parent] = y;
            } else {
                right_child[x_parent] = y;
            }
        } else {
            root = y;
            parent.erase(y);
        }

        // Set y as parent of x
        left_child[y] = x;
        parent[x] = y;

        // fix ranges up from y
        updateRangesUp(y);
    }

    // Rotate x and its left child up (right-rotate)
    void rotateRight(int x) {
        if (x == -1) return;

        auto left_it = left_child.find(x);
        if (left_it == left_child.end()) return; // no left child

        int y = left_it->second;
        if (dummy_leaves.find(y) != dummy_leaves.end()) return; // can't rotate dummy

        // Perform rotation logic
        auto y_right_it = right_child.find(y);
        int y_right = (y_right_it != right_child.end()) ? y_right_it->second : -1;

        // Update x's left child
        if (y_right != -1) {
            left_child[x] = y_right;
            parent[y_right] = x;
        } else {
            left_child.erase(x);
        }

        // Update y's parent
        auto x_parent_it = parent.find(x);
        if (x_parent_it != parent.end()) {
            int x_parent = x_parent_it->second;
            parent[y] = x_parent;

            auto x_parent_left_it = left_child.find(x_parent);
            if (x_parent_left_it != left_child.end() && x_parent_left_it->second == x) {
                left_child[x_parent] = y;
            } else {
                right_child[x_parent] = y;
            }
        } else {
            root = y;
            parent.erase(y);
        }

        // Set y as parent of x
        right_child[y] = x;
        parent[x] = y;

        // fix ranges up from y
        updateRangesUp(y);
    }

    // Collect all parent->child edges into `out_set`
    void collectEdges(int n,
                      std::unordered_set<std::pair<int,int>, PairHash, PairEq>& out_set) const
    {
        if (n == -1) return;

        auto left_it = left_child.find(n);
        if (left_it != left_child.end() && dummy_leaves.find(left_it->second) == dummy_leaves.end()) {
            out_set.insert({n, left_it->second});
            collectEdges(left_it->second, out_set);
        }

        auto right_it = right_child.find(n);
        if (right_it != right_child.end() && dummy_leaves.find(right_it->second) == dummy_leaves.end()) {
            out_set.insert({n, right_it->second});
            collectEdges(right_it->second, out_set);
        }
    }

    // Print ranges for all original nodes
    void print() const {
        std::cout << "Tree structure:\n";
        printTreeStructure();

        std::cout << "\nRanges (original nodes only):\n";

        // Get all original nodes and sort them
        std::vector<int> nodes_vec;
        for (const auto& [node, range] : ranges) {
            if (dummy_leaves.find(node) == dummy_leaves.end()) {
                nodes_vec.push_back(node);
            }
        }
        std::sort(nodes_vec.begin(), nodes_vec.end());

        // Print all original nodes
        for (int node : nodes_vec) {
            auto range = ranges.at(node);
            std::cout << node << "(" << range.first << "," << range.second << ") ";
        }
        std::cout << std::endl;
    }

    // Helper functions for tree navigation
    int getParent(int node) const {
        auto it = parent.find(node);
        return (it != parent.end()) ? it->second : -1;
    }

    int getLeftChild(int node) const {
        auto it = left_child.find(node);
        return (it != left_child.end()) ? it->second : -1;
    }

    int getRightChild(int node) const {
        auto it = right_child.find(node);
        return (it != right_child.end()) ? it->second : -1;
    }

private:
    // Build leaves vector in left-to-right order
    void buildLeavesVector(int node) {
        if (node == -1) return;

        auto left_it = left_child.find(node);
        auto right_it = right_child.find(node);

        // If it's a leaf (no children), add to leaves
        if (left_it == left_child.end() && right_it == right_child.end()) {
            leaves.push_back(node);
            return;
        }

        // Recursively build leaves in inorder
        if (left_it != left_child.end()) {
            buildLeavesVector(left_it->second);
        }
        if (right_it != right_child.end()) {
            buildLeavesVector(right_it->second);
        }
    }

    // Build tree recursively
    void buildRecursive(const std::vector<int>& preorder, const std::vector<int>& inorder,
                        int pre_start, int pre_end, int in_start, int in_end, int par) {
        if (pre_start > pre_end || in_start > in_end) {
            return;
        }

        int root_val = preorder[pre_start];

        // Set parent relationship
        if (par != -1) {
            parent[root_val] = par;
        }

        // Find root position in inorder
        int root_idx = -1;
        for (int i = in_start; i <= in_end; i++) {
            if (inorder[i] == root_val) {
                root_idx = i;
                break;
            }
        }

        int left_size = root_idx - in_start;

        // Build left subtree
        if (left_size > 0) {
            int left_root = preorder[pre_start + 1];
            left_child[root_val] = left_root;
            buildRecursive(preorder, inorder,
                           pre_start + 1, pre_start + left_size,
                           in_start, root_idx - 1, root_val);
        }

        // Build right subtree
        if (root_idx < in_end) {
            int right_root = preorder[pre_start + left_size + 1];
            right_child[root_val] = right_root;
            buildRecursive(preorder, inorder,
                           pre_start + left_size + 1, pre_end,
                           root_idx + 1, in_end, root_val);
        }
    }

    // Add dummy leaves to make every original node an inner node
    void addDummyLeaves() {
        std::vector<int> original_nodes;

        // Collect all original nodes
        for (int node : original_inorder) {
            original_nodes.push_back(node);
        }

        // Add dummy leaves to nodes that don't have both children
        for (int node : original_nodes) {
            int left_child_node = getLeftChild(node);
            int right_child_node = getRightChild(node);

            // Add dummy left child if missing
            if (left_child_node == -1) {
                int dummy_left = next_dummy_id--;
                dummy_leaves.insert(dummy_left);
                left_child[node] = dummy_left;
                parent[dummy_left] = node;
            }

            // Add dummy right child if missing
            if (right_child_node == -1) {
                int dummy_right = next_dummy_id--;
                dummy_leaves.insert(dummy_right);
                right_child[node] = dummy_right;
                parent[dummy_right] = node;
            }
        }

        // Recursively add dummy leaves to newly created inner nodes if needed
        addDummyLeavesRecursive();
    }

    // Recursively ensure all non-dummy nodes have two children
    void addDummyLeavesRecursive() {
        bool added_new_dummies = true;

        while (added_new_dummies) {
            added_new_dummies = false;
            std::vector<int> nodes_to_check;

            // Collect all non-dummy nodes
            for (const auto& [child, par] : parent) {
                if (dummy_leaves.find(child) == dummy_leaves.end()) {
                    nodes_to_check.push_back(child);
                }
                if (dummy_leaves.find(par) == dummy_leaves.end()) {
                    nodes_to_check.push_back(par);
                }
            }

            // Add root if it exists
            if (root != -1 && dummy_leaves.find(root) == dummy_leaves.end()) {
                nodes_to_check.push_back(root);
            }

            // Remove duplicates
            std::sort(nodes_to_check.begin(), nodes_to_check.end());
            nodes_to_check.erase(std::unique(nodes_to_check.begin(), nodes_to_check.end()), nodes_to_check.end());

            for (int node : nodes_to_check) {
                if (dummy_leaves.find(node) != dummy_leaves.end()) continue;

                int left_child_node = getLeftChild(node);
                int right_child_node = getRightChild(node);

                // Add dummy left child if missing
                if (left_child_node == -1) {
                    int dummy_left = next_dummy_id--;
                    dummy_leaves.insert(dummy_left);
                    left_child[node] = dummy_left;
                    parent[dummy_left] = node;
                    added_new_dummies = true;
                }

                // Add dummy right child if missing
                if (right_child_node == -1) {
                    int dummy_right = next_dummy_id--;
                    dummy_leaves.insert(dummy_right);
                    right_child[node] = dummy_right;
                    parent[dummy_right] = node;
                    added_new_dummies = true;
                }
            }
        }
    }

    // TRUE O(1) range update for a single node (with dummy leaves)
    void updateNodeRange(int node) {
        if (dummy_leaves.find(node) != dummy_leaves.end()) {
            // Dummy leaves have empty ranges
            ranges[node] = {0, 0};
            return;
        }

        // For inner nodes with dummy leaves, range calculation is O(1)
        int left_child_node = getLeftChild(node);
        int right_child_node = getRightChild(node);

        // With dummy leaves, every original node is guaranteed to have both children
        if (left_child_node == -1 || right_child_node == -1) {
            // This shouldn't happen with proper dummy leaf augmentation
            ranges[node] = {position_map[node], position_map[node] + 1};
            return;
        }

        // O(1) calculation: combine ranges from left and right children
        auto left_range = ranges[left_child_node];
        auto right_range = ranges[right_child_node];

        int start = std::min(left_range.first, position_map[node]);
        int end = std::max(right_range.second, position_map[node] + 1);

        // If left child is dummy, start from this node's position
        if (dummy_leaves.find(left_child_node) != dummy_leaves.end()) {
            start = position_map[node];
        }

        // If right child is dummy, end at this node's position + 1
        if (dummy_leaves.find(right_child_node) != dummy_leaves.end()) {
            end = position_map[node] + 1;
        }

        ranges[node] = {start, end};
    }

    // Walk up from n to root, fixing ranges
    void updateRangesUp(int n) {
        while (n != -1) {
            updateNodeRange(n);
            auto parent_it = parent.find(n);
            n = (parent_it != parent.end()) ? parent_it->second : -1;
        }
    }

    // Calculate ranges for all nodes - only called during initial build
    void calculateAllRanges() {
        ranges.clear();

        // Initialize dummy leaf ranges (empty ranges)
        for (int dummy : dummy_leaves) {
            ranges[dummy] = {0, 0};
        }

        // Calculate ranges for original nodes in post-order (bottom-up)
        calculateRangesPostOrder(root);
    }

    // Post-order traversal to calculate ranges bottom-up
    void calculateRangesPostOrder(int node) {
        if (node == -1) return;

        int left_child_node = getLeftChild(node);
        int right_child_node = getRightChild(node);

        // Process children first (post-order)
        if (left_child_node != -1) {
            calculateRangesPostOrder(left_child_node);
        }
        if (right_child_node != -1) {
            calculateRangesPostOrder(right_child_node);
        }

        // Now calculate this node's range
        updateNodeRange(node);
    }
};

// Compute and print shared edges between two RangeTreeMap trees
void getSharedEdges(const RangeTreeMap& A, const RangeTreeMap& B) {
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
    RangeTreeMap tree;

    struct TestCase { std::string name;
        std::vector<int> preorder;
        std::vector<int> inorder;
        int rotateAt;
        bool doLeft; };

    std::vector<TestCase> tests = {
            {"Test 1: Simple 3-node, left-rotate at 2", {2,1,3}, {1,2,3}, 2, true},
            {"Test 2: Right-chain of 4, left-rotate at 2", {1,2,3,4}, {1,2,3,4}, 2, true},
            {"Test 3: Balanced 7-node, right-rotate at 6", {4,2,1,3,6,5,7}, {1,2,3,4,5,6,7}, 6, false}
    };

    for(auto& tc: tests) {
        std::cout << tc.name << "\n";
        tree.build(tc.preorder, tc.inorder);
        std::cout << "Before:\n";
        tree.printTreeStructure();
        tree.printLeaves();
        std::cout << "\n";
        int x = tc.rotateAt; // In map version, we use the value directly
        if(tc.doLeft) tree.rotateLeft(x); else tree.rotateRight(x);
        std::cout << "After:\n";
        tree.printTreeStructure();
        tree.printLeaves();
        std::cout << "\n\n";
    }

    // Test 4: Multiple rotations on 7-node
    std::cout << "Test 4: Multiple rotations on 7-node\n";
    tree.build({4,2,1,3,6,5,7}, {1,2,3,4,5,6,7});
    std::cout << "Original structure:\n";
    tree.printTreeStructure();
    tree.printLeaves();
    std::cout << "\n";
    tree.rotateLeft(2);
    tree.rotateRight(6);
    std::cout << "After rotateLeft(2) and rotateRight(6):\n";
    tree.printTreeStructure();
    tree.printLeaves();
    std::cout << "\n";

    // Test 5: rotateLeft on node with no right child
    std::cout << "Test 5: rotateLeft on leaf node (1)\n";
    tree.build({2,1,3}, {1,2,3});
    std::cout << "Before (should be unchanged):\n";
    tree.printTreeStructure();
    tree.printLeaves();
    std::cout << "\n";
    tree.rotateLeft(1); // 1 has no right child
    std::cout << "After (should still be unchanged):\n";
    tree.printTreeStructure();
    tree.printLeaves();
    std::cout << "\n";

    // Test 6: Shared-edge detection between two trees
    std::cout << "Test 6: Shared edges between two trees\n";
    // Build a second tree to compare
    RangeTreeMap tree2;
    tree.build({4,2,6,5,7}, {2,4,5,6,7}); // Tree1
    tree2.build({4,2,6,3,5}, {2,4,3,6,5}); // Tree2
    std::cout << "Tree1 structure:\n";
    tree.printTreeStructure();
    tree.printLeaves();
    std::cout << "\nTree2 structure:\n";
    tree2.printTreeStructure();
    tree2.printLeaves();
    std::cout << "\nShared edges:\n";
    getSharedEdges(tree, tree2);

    std::cout << "\n Additional Tests for shared edge\n\n";

    // Both trees empty
    {
        std::cout << "-- Case 1: Both trees empty --\n";
        RangeTreeMap A, B;
        A.build({}, {});
        B.build({}, {});
        std::cout << "Tree A and B are both empty.\n";
        std::cout << "Shared edges (expected none):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // One tree empty, one nonempty
    {
        std::cout << "Case 2: One empty, one nonempty\n";
        RangeTreeMap A, B;
        A.build({}, {});                    // A is empty
        B.build({1,2,3}, {2,1,3});          // B:   1
        //      / \
                                            //     2   3
        std::cout << "Tree A (empty), Tree B:\n";
        std::cout << "Tree B structure:\n";
        B.printTreeStructure();
        B.printLeaves();
        std::cout << "Shared edges (expected none):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // Identical small tree
    {
        std::cout << "Case 3: Identical small tree\n";
        RangeTreeMap A, B;
        A.build({2,1,3}, {1,2,3});  //      2
        //     / \
                                    //    1   3
        B.build({2,1,3}, {1,2,3});
        std::cout << "Tree A and Tree B (identical):\n";
        A.printTreeStructure();
        A.printLeaves();
        std::cout << "Shared edges (expected: 2->1, 2->3):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // Same values, different shape
    {
        std::cout << "Case 4: Same values, different shape\n";
        RangeTreeMap A, B;
        A.build({4,2,3,1}, {2,3,1,4});
        // Tree A:
        //    4
        //   /
        //  2
        //   \
        //    3
        //     \
        //      1
        std::cout << "Tree A structure:\n";
        A.printTreeStructure();
        A.printLeaves();

        B.build({4,3,2,1}, {1,2,3,4});
        // Tree B:
        //    4
        //     \
        //      3
        //     /
        //    2
        //   /
        //  1
        std::cout << "Tree B structure:\n";
        B.printTreeStructure();
        B.printLeaves();

        std::cout << "Shared edges (expected none):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // Partial overlap deep inside
    {
        std::cout << "Case 5: Partial overlap deep inside\n";
        RangeTreeMap A, B;
        A.build({5,3,1,8,7,9}, {1,3,5,7,8,9});
        // Tree A:
        //      5
        //     / \
        //    3   8
        //   /   / \
        //  1   7   9
        std::cout << "Tree A structure:\n";
        A.printTreeStructure();
        A.printLeaves();

        B.build({5,2,8,7}, {2,5,7,8});
        // Tree B:
        //      5
        //     / \
        //    2   8
        //       /
        //      7
        std::cout << "Tree B structure:\n";
        B.printTreeStructure();
        B.printLeaves();

        std::cout << "Shared edges (expected: 5->8, 8->7):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    return 0;
}