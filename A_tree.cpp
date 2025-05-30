#include <iostream>
#include <vector>
#include <unordered_map>
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

    RangeTree() : root(-1) {}

    // Build tree from preorder and inorder traversals
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        ranges.clear();
        parent.clear();
        left_child.clear();
        right_child.clear();
        root = -1;

        // Store the original inorder sequence
        original_inorder = inorder;

        if (!preorder.empty()) {
            root = preorder[0];
            buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, -1);
            calculateAllRanges();
        }
    }

    void print() const {
        std::cout << "Tree structure:\n";
        printStructure();

        std::cout << "\nRanges (actual nodes only):\n";

        // Get all nodes and sort them
        std::vector<int> nodes;
        for (const auto& [node, range] : ranges) {
            nodes.push_back(node);
        }
        std::sort(nodes.begin(), nodes.end());

        // Print all nodes
        for (int node : nodes) {
            auto range = ranges.at(node);
            std::cout << node << "(" << range.first << "," << range.second << ") ";
        }
        std::cout << std::endl;
    }

    // Print tree structure for debugging
    void printStructure() const {
        std::cout << "Root: " << root << "\n";
        std::cout << "Parent relationships:\n";
        for (const auto& [child, par] : parent) {
            std::cout << "  " << child << " -> parent: " << par << "\n";
        }
        std::cout << "Child relationships:\n";
        for (const auto& [par, left] : left_child) {
            int right = (right_child.find(par) != right_child.end()) ? right_child.at(par) : -1;
            std::cout << "  " << par << " -> left: " << left << ", right: " << right << "\n";
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

    // Perform right rotation at the given node
    bool rotateRight(int node) {
        int left = getLeftChild(node);
        if (left == -1) return false;

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

        // Recalculate ranges after rotation
        calculateAllRanges();
        return true;
    }

    // Perform left rotation at the given node
    bool rotateLeft(int node) {
        int right = getRightChild(node);
        if (right == -1) return false;

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

        // Recalculate ranges after rotation
        calculateAllRanges();
        return true;
    }

private:
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

    // Calculate ranges for all nodes based on ORIGINAL inorder sequence
    void calculateAllRanges() {
        ranges.clear();

        // Create position mapping from ORIGINAL inorder sequence
        std::unordered_map<int, int> position_map;
        for (size_t i = 0; i < original_inorder.size(); i++) {
            position_map[original_inorder[i]] = i;
        }

        // Calculate range for each node based on its current subtree span
        for (int node : original_inorder) {
            auto [leftmost_pos, rightmost_pos] = getSubtreePositionRange(node, position_map);
            ranges[node] = {leftmost_pos, rightmost_pos + 1}; // +1 for exclusive end
        }
    }

    // Get the leftmost and rightmost positions that a subtree spans
    std::pair<int, int> getSubtreePositionRange(int node, const std::unordered_map<int, int>& position_map) {
        if (node == -1) return {-1, -1};

        // Start with this node's position in the ORIGINAL sequence
        int node_pos = position_map.at(node);
        int leftmost = node_pos;
        int rightmost = node_pos;

        // Check left subtree
        int left_child_node = getLeftChild(node);
        if (left_child_node != -1) {
            auto [left_min, left_max] = getSubtreePositionRange(left_child_node, position_map);
            leftmost = std::min(leftmost, left_min);
            rightmost = std::max(rightmost, left_max);
        }

        // Check right subtree
        int right_child_node = getRightChild(node);
        if (right_child_node != -1) {
            auto [right_min, right_max] = getSubtreePositionRange(right_child_node, position_map);
            leftmost = std::min(leftmost, right_min);
            rightmost = std::max(rightmost, right_max);
        }

        return {leftmost, rightmost};
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
    std::cout << "=== TESTING BASIC RIGHT ROTATION ===\n";

    RangeTree tree = buildExampleTree();
    std::cout << "Original tree:\n";
    tree.print();

    std::cout << "\nAfter right rotation at node 4:\n";
    tree.rotateRight(4);
    tree.print();

    std::cout << "\nExpected ranges:\n";
    std::cout << "1(0,1) 2(0,2) 3(0,5) 4(3,5) 5(4,5) 6(0,6)\n";
}

void testLeftRotation() {
    std::cout << "\n=== TESTING LEFT ROTATION ===\n";

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
    std::cout << "\n=== TESTING MULTIPLE ROTATIONS ===\n";

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