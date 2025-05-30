#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>

struct RangeTree {
    std::unordered_map<int, std::pair<int, int>> ranges; // node_value -> (start, end) range

    // Parent-child relationships
    std::unordered_map<int, int> parent;       // child -> parent
    std::unordered_map<int, int> left_child;   // parent -> left child
    std::unordered_map<int, int> right_child;  // parent -> right child
    int root;                                  // root node value

    // Store the original inorder sequence for range calculations
    std::vector<int> original_inorder;

    // Set to track which nodes are dummy leaves (negative values)
    std::unordered_set<int> dummy_leaves;

    // Position mapping for O(1) range calculations
    std::unordered_map<int, int> position_map;

    // Next dummy leaf ID (negative numbers)
    int next_dummy_id;

    RangeTree() : root(-1), next_dummy_id(-1000) {}

    // Build tree from preorder and inorder traversals
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        ranges.clear();
        parent.clear();
        left_child.clear();
        right_child.clear();
        dummy_leaves.clear();
        position_map.clear();
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
            addDummyLeaves();
            calculateAllRanges();
        }
    }

    void print() const {
        std::cout << "Tree structure:\n";
        printStructure();

        std::cout << "\nRanges (original nodes only):\n";

        // Get all original nodes and sort them
        std::vector<int> nodes;
        for (const auto& [node, range] : ranges) {
            if (dummy_leaves.find(node) == dummy_leaves.end()) {
                nodes.push_back(node);
            }
        }
        std::sort(nodes.begin(), nodes.end());

        // Print all original nodes
        for (int node : nodes) {
            auto range = ranges.at(node);
            std::cout << node << "(" << range.first << "," << range.second << ") ";
        }
        std::cout << std::endl;
    }

    // Print tree structure for debugging
    void printStructure() const {
        std::cout << "Root: " << root << "\n";
        std::cout << "Parent relationships (including dummy leaves):\n";
        for (const auto& [child, par] : parent) {
            std::string child_type = (dummy_leaves.find(child) != dummy_leaves.end()) ? " (dummy)" : "";
            std::cout << "  " << child << child_type << " -> parent: " << par << "\n";
        }
        std::cout << "Child relationships:\n";
        for (const auto& [par, left] : left_child) {
            int right = (right_child.find(par) != right_child.end()) ? right_child.at(par) : -1;
            std::string left_type = (left != -1 && dummy_leaves.find(left) != dummy_leaves.end()) ? " (dummy)" : "";
            std::string right_type = (right != -1 && dummy_leaves.find(right) != dummy_leaves.end()) ? " (dummy)" : "";
            std::cout << "  " << par << " -> left: " << left << left_type << ", right: " << right << right_type << "\n";
        }
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

    // Perform right rotation at the given node - TRUE O(1)!
    bool rotateRight(int node) {
        int left = getLeftChild(node);
        if (left == -1 || dummy_leaves.find(left) != dummy_leaves.end()) return false;

        int par = getParent(node);
        int left_right = getRightChild(left);

        // Update parent relationships
        if (par != -1) {
            if (getLeftChild(par) == node) {
                left_child[par] = left;
            } else {
                right_child[par] = left;
            }
            parent[left] = par;
        } else {
            root = left;
            parent.erase(left);
        }

        parent[node] = left;
        right_child[left] = node;

        if (left_right != -1) {
            parent[left_right] = node;
            left_child[node] = left_right;
        } else {
            left_child.erase(node);
        }

        // TRUE O(1) range updates - only update the two affected nodes!
        updateNodeRange(node);  // The node that was rotated down
        updateNodeRange(left);  // The node that was rotated up

        return true;
    }

    // Perform left rotation at the given node - TRUE O(1)!
    bool rotateLeft(int node) {
        int right = getRightChild(node);
        if (right == -1 || dummy_leaves.find(right) != dummy_leaves.end()) return false;

        int par = getParent(node);
        int right_left = getLeftChild(right);

        // Update parent relationships
        if (par != -1) {
            if (getLeftChild(par) == node) {
                left_child[par] = right;
            } else {
                right_child[par] = right;
            }
            parent[right] = par;
        } else {
            root = right;
            parent.erase(right);
        }

        parent[node] = right;
        left_child[right] = node;

        if (right_left != -1) {
            parent[right_left] = node;
            right_child[node] = right_left;
        } else {
            right_child.erase(node);
        }

        // TRUE O(1) range updates - only update the two affected nodes!
        updateNodeRange(node);   // The node that was rotated down
        updateNodeRange(right);  // The node that was rotated up

        return true;
    }

private:
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
        // (This handles cases where we need to add more levels)
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

// Build the example tree from your professor
RangeTree buildExampleTree() {
    RangeTree tree;

    // Building the example tree:
    //         6
    //      /
    //     4
    //     /\
    //   3  5
    //  /
    // 2
    // /
    // 1

    // Preorder: [6, 4, 3, 2, 1, 5]
    // Inorder: [1, 2, 3, 4, 5, 6]

    std::vector<int> preorder = {6, 4, 3, 2, 1, 5};
    std::vector<int> inorder = {1, 2, 3, 4, 5, 6};

    tree.build(preorder, inorder);

    return tree;
}

// Build a balanced binary tree for testing
RangeTree buildBalancedTree() {
    RangeTree tree;

    // Building a balanced tree:
    //       4
    //     /   \
    //    2     6
    //   / \   / \
    //  1   3 5   7

    // Preorder: [4, 2, 1, 3, 6, 5, 7]
    // Inorder: [1, 2, 3, 4, 5, 6, 7]

    std::vector<int> preorder = {4, 2, 1, 3, 6, 5, 7};
    std::vector<int> inorder = {1, 2, 3, 4, 5, 6, 7};

    tree.build(preorder, inorder);

    return tree;
}

void testBasicRotation() {
    std::cout << "=== TESTING BASIC RIGHT ROTATION WITH DUMMY LEAVES ===\n";

    RangeTree tree = buildExampleTree();
    std::cout << "Original tree (with dummy leaves added):\n";
    tree.print();

    std::cout << "\nAfter right rotation at node 4:\n";
    tree.rotateRight(4);
    tree.print();

    std::cout << "\nExpected ranges:\n";
    std::cout << "1(0,1) 2(0,2) 3(0,5) 4(3,5) 5(4,5) 6(0,6)\n";
}

void testLeftRotation() {
    std::cout << "\n=== TESTING LEFT ROTATION WITH DUMMY LEAVES ===\n";

    // Create a tree where left rotation is possible
    RangeTree tree;
    std::vector<int> preorder = {2, 1, 4, 3, 5};
    std::vector<int> inorder = {1, 2, 3, 4, 5};
    tree.build(preorder, inorder);

    std::cout << "Tree before left rotation at node 2:\n";
    tree.print();

    std::cout << "\nAfter left rotation at node 2:\n";
    tree.rotateLeft(2);
    tree.print();
}

void testMultipleRotations() {
    std::cout << "\n=== TESTING MULTIPLE ROTATIONS WITH DUMMY LEAVES ===\n";

    RangeTree tree = buildBalancedTree();
    std::cout << "Balanced tree:\n";
    tree.print();

    std::cout << "\nAfter right rotation at root (4):\n";
    tree.rotateRight(4);
    tree.print();

    std::cout << "\nAfter left rotation at new root (2):\n";
    tree.rotateLeft(2);
    tree.print();
    std::cout << "Should restore original structure\n";
}

int main() {
    testBasicRotation();
    testLeftRotation();
    testMultipleRotations();
    return 0;
}