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

    RangeTree() : root(-1) {}

    // Build tree from preorder and inorder traversals
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        ranges.clear();
        parent.clear();
        left_child.clear();
        right_child.clear();
        root = -1;

        if (!preorder.empty()) {
            root = preorder[0];
            buildRecursive(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, -1);
            calculateAllRanges(inorder);
        }
    }

    // Print the tree representation
    void print() const {
        std::cout << "Tree structure:\n";
        printStructure();

        std::cout << "\nRanges (actual nodes only):\n";
        std::vector<int> node_order = {1, 2, 3, 4, 5, 6}; // Adjust based on your tree
        for (int node : node_order) {
            if (ranges.find(node) != ranges.end()) {
                auto range = ranges.at(node);
                std::cout << node << "(" << range.first << "," << range.second << ") ";
            }
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
        if (left == -1) return false; // Can't rotate without left child

        int par = getParent(node);
        int left_right = getRightChild(left);

        // Update parent relationships
        if (par != -1) {
            // Update parent's child pointer
            if (getLeftChild(par) == node) {
                left_child[par] = left;
            } else {
                right_child[par] = left;
            }
            parent[left] = par;
        } else {
            // Node was root, now left is root
            root = left;
            parent.erase(left);
        }

        // Left becomes parent of node
        parent[node] = left;
        right_child[left] = node;

        // Left's right child becomes node's left child
        if (left_right != -1) {
            parent[left_right] = node;
            left_child[node] = left_right;
        } else {
            left_child.erase(node);
        }

        // Recalculate ranges for the entire tree
        std::vector<int> inorder_seq;
        inorderTraversal(root, inorder_seq);
        calculateAllRanges(inorder_seq);

        return true;
    }

    // Perform left rotation at the given node
    bool rotateLeft(int node) {
        int right = getRightChild(node);
        if (right == -1) return false; // Can't rotate without right child

        int par = getParent(node);
        int right_left = getLeftChild(right);

        // Update parent relationships
        if (par != -1) {
            // Update parent's child pointer
            if (getLeftChild(par) == node) {
                left_child[par] = right;
            } else {
                right_child[par] = right;
            }
            parent[right] = par;
        } else {
            // Node was root, now right is root
            root = right;
            parent.erase(right);
        }

        // Right becomes parent of node
        parent[node] = right;
        left_child[right] = node;

        // Right's left child becomes node's right child
        if (right_left != -1) {
            parent[right_left] = node;
            right_child[node] = right_left;
        } else {
            right_child.erase(node);
        }

        // Recalculate ranges for the entire tree
        std::vector<int> inorder_seq;
        inorderTraversal(root, inorder_seq);
        calculateAllRanges(inorder_seq);

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

    // Calculate ranges for all nodes based on inorder traversal
    void calculateAllRanges(const std::vector<int>& inorder_seq) {
        ranges.clear();

        // For each node, find its position and subtree coverage in inorder
        for (int node : inorder_seq) {
            // Find the leftmost and rightmost positions this node covers
            int left_pos = findLeftmostInSubtree(node, inorder_seq);
            int right_pos = findRightmostInSubtree(node, inorder_seq);
            ranges[node] = {left_pos, right_pos + 1}; // +1 for exclusive end
        }
    }

    // Find leftmost position in subtree rooted at node
    int findLeftmostInSubtree(int node, const std::vector<int>& inorder_seq) {
        // Get all nodes in subtree
        std::vector<int> subtree_nodes;
        getSubtreeInorder(node, subtree_nodes);

        // Find the position of the leftmost node in the full inorder sequence
        if (subtree_nodes.empty()) return 0;

        int leftmost_node = subtree_nodes[0];
        for (size_t i = 0; i < inorder_seq.size(); i++) {
            if (inorder_seq[i] == leftmost_node) {
                return i;
            }
        }
        return 0;
    }

    // Find rightmost position in subtree rooted at node
    int findRightmostInSubtree(int node, const std::vector<int>& inorder_seq) {
        // Get all nodes in subtree
        std::vector<int> subtree_nodes;
        getSubtreeInorder(node, subtree_nodes);

        // Find the position of the rightmost node in the full inorder sequence
        if (subtree_nodes.empty()) return 0;

        int rightmost_node = subtree_nodes.back();
        for (size_t i = 0; i < inorder_seq.size(); i++) {
            if (inorder_seq[i] == rightmost_node) {
                return i;
            }
        }
        return 0;
    }

    // Get inorder traversal of subtree rooted at node
    void getSubtreeInorder(int node, std::vector<int>& result) {
        if (node == -1) return;

        getSubtreeInorder(getLeftChild(node), result);
        result.push_back(node);
        getSubtreeInorder(getRightChild(node), result);
    }

    // Inorder traversal to get nodes in left-to-right order
    void inorderTraversal(int node, std::vector<int>& result) {
        if (node == -1) return;

        inorderTraversal(getLeftChild(node), result);
        result.push_back(node);
        inorderTraversal(getRightChild(node), result);
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

int main() {
    std::cout << "=== ORIGINAL TREE ===\n";
    RangeTree tree = buildExampleTree();
    tree.print();

    std::cout << "\n=== TESTING RIGHT ROTATION AT NODE 4 ===\n";

    // Make a copy for rotation
    RangeTree rotated_tree = tree;

    std::cout << "Before rotation:\n";
    rotated_tree.print();

    // Perform right rotation at node 4
    if (rotated_tree.rotateRight(4)) {
        std::cout << "\nAfter right rotation at node 4:\n";
        rotated_tree.print();

        std::cout << "\nExpected after rotation:\n";
        std::cout << "Structure: 6->3, 3->2,4, 2->1, 4->5\n";
        std::cout << "Expected ranges: 1(0,1) 2(0,2) 3(0,5) 4(3,5) 5(4,5) 6(0,6)\n";

    } else {
        std::cout << "Right rotation at node 4 failed!\n";
    }

    std::cout << "\n=== TESTING LEFT ROTATION AT NODE 4 ===\n";

    // Test left rotation on original tree
    RangeTree left_rotated_tree = tree;

    if (left_rotated_tree.rotateLeft(4)) {
        std::cout << "After left rotation at node 4:\n";
        left_rotated_tree.print();
    } else {
        std::cout << "Left rotation at node 4 failed (expected - no right child initially)!\n";
    }

    return 0;
}