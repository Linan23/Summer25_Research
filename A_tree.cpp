#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>

struct RangeTree {
    std::vector<int> leaves;                   // Store leaf values in left-to-right order
    std::unordered_map<int, std::pair<int, int>> ranges; // node_value -> (start, end) leaf range

    // Parent-child relationships
    std::unordered_map<int, int> parent;       // child -> parent
    std::unordered_map<int, int> left_child;   // parent -> left child
    std::unordered_map<int, int> right_child;  // parent -> right child
    int root;                                  // root node value

    RangeTree() : root(-1) {}

    // Build tree and calculate ranges based on leaf positions
    void buildWithLeaves(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        leaves.clear();
        ranges.clear();
        parent.clear();
        left_child.clear();
        right_child.clear();
        root = -1;

        if (!preorder.empty()) {
            root = preorder[0];
            buildRecursiveWithLeaves(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1, -1);
        }
    }

    // Print the tree representation
    void print() const {
        std::cout << "Tree structure:\n";
        printStructure();

        std::cout << "\nLeaves in left-to-right order: ";
        for (size_t i = 0; i < leaves.size(); i++) {
            std::cout << "L" << i << " ";
        }

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

    // Check if two trees have the same structure (for comparison)
    bool hasSameStructure(const RangeTree& other) const {
        // Compare ranges - if ranges match, trees are equivalent
        if (ranges.size() != other.ranges.size()) return false;

        for (const auto& [node, range] : ranges) {
            auto it = other.ranges.find(node);
            if (it == other.ranges.end() || it->second != range) {
                return false;
            }
        }
        return true;
    }

    // Get the total number of leaves
    int leafCount() const {
        return leaves.size();
    }

    // Perform right rotation at the given node
    // Right rotation: if node has left child, promote left child up
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

        // Recalculate ranges for affected subtrees
        recalculateRanges();

        return true;
    }

    // Perform left rotation at the given node
    // Left rotation: if node has right child, promote right child up
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

        // Recalculate ranges for affected subtrees
        recalculateRanges();

        return true;
    }

    // Get all possible rotations from current state (for brute force)
    std::vector<RangeTree> getAllPossibleRotations() const {
        std::vector<RangeTree> rotations;

        // Try rotating at every node
        for (const auto& [node, range] : ranges) {
            // Try right rotation
            RangeTree right_rotated = *this;
            if (right_rotated.rotateRight(node)) {
                rotations.push_back(right_rotated);
            }

            // Try left rotation
            RangeTree left_rotated = *this;
            if (left_rotated.rotateLeft(node)) {
                rotations.push_back(left_rotated);
            }
        }

        return rotations;
    }

private:
    // Recalculate ranges for all nodes after a rotation
    void recalculateRanges() {
        ranges.clear();
        if (root != -1) {
            calculateRangesRecursive(root);
        }
    }

    // Recursively calculate ranges for each node based on current tree structure
    std::pair<int, int> calculateRangesRecursive(int node) {
        int left = getLeftChild(node);
        int right = getRightChild(node);

        std::pair<int, int> left_range, right_range;

        // Calculate left subtree range
        if (left != -1) {
            left_range = calculateRangesRecursive(left);
        } else {
            // No left child - get the leaf position that would be here
            int leaf_pos = getLeafPosition(node, true);
            left_range = {leaf_pos, leaf_pos};
        }

        // Calculate right subtree range
        if (right != -1) {
            right_range = calculateRangesRecursive(right);
        } else {
            // No right child - get the leaf position that would be here
            int leaf_pos = getLeafPosition(node, false);
            right_range = {leaf_pos, leaf_pos};
        }

        // This node's range spans from leftmost to rightmost leaf in its subtree
        int range_start = left_range.first;
        int range_end = right_range.second;

        ranges[node] = {range_start, range_end};
        return {range_start, range_end};
    }

    // Get the leaf position for a missing child
    int getLeafPosition(int node, bool isLeft) {
        // Perform inorder traversal to get the sequence of nodes
        std::vector<int> inorder_nodes;
        inorderTraversal(root, inorder_nodes);

        // Find position of current node
        int node_index = -1;
        for (size_t i = 0; i < inorder_nodes.size(); i++) {
            if (inorder_nodes[i] == node) {
                node_index = i;
                break;
            }
        }

        if (isLeft) {
            // Left missing child corresponds to leaf position before this node
            return node_index;
        } else {
            // Right missing child corresponds to leaf position after this node
            return node_index + 1;
        }
    }

    // Inorder traversal to get nodes in left-to-right order
    void inorderTraversal(int node, std::vector<int>& result) {
        if (node == -1) return;

        int left = getLeftChild(node);
        int right = getRightChild(node);

        if (left != -1) {
            inorderTraversal(left, result);
        }

        result.push_back(node);

        if (right != -1) {
            inorderTraversal(right, result);
        }
    }

    std::pair<int, int> buildRecursiveWithLeaves(const std::vector<int>& preorder, const std::vector<int>& inorder,
                                                 int pre_start, int pre_end, int in_start, int in_end, int par) {
        if (pre_start > pre_end || in_start > in_end) {
            // Add a leaf and return its position
            int leaf_pos = leaves.size();
            leaves.push_back(1000 + leaf_pos); // Unique leaf value
            return {leaf_pos, leaf_pos};
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
        int right_size = in_end - root_idx;

        // Build left subtree
        std::pair<int, int> left_range;
        if (left_size > 0) {
            int left_root = preorder[pre_start + 1];
            left_child[root_val] = left_root;
            left_range = buildRecursiveWithLeaves(preorder, inorder,
                                                  pre_start + 1, pre_start + left_size,
                                                  in_start, root_idx - 1, root_val);
        } else {
            // No left child - range represents a leaf
            left_range = buildRecursiveWithLeaves(preorder, inorder,
                                                  pre_start + 1, pre_start,
                                                  in_start, root_idx - 1, root_val);
        }

        // Build right subtree
        std::pair<int, int> right_range;
        if (right_size > 0) {
            int right_root = preorder[pre_start + left_size + 1];
            right_child[root_val] = right_root;
            right_range = buildRecursiveWithLeaves(preorder, inorder,
                                                   pre_start + left_size + 1, pre_end,
                                                   root_idx + 1, in_end, root_val);
        } else {
            // No right child - range represents a leaf
            right_range = buildRecursiveWithLeaves(preorder, inorder,
                                                   pre_start + left_size + 1, pre_start + left_size,
                                                   root_idx + 1, in_end, root_val);
        }

        // This node's range spans from leftmost leaf to rightmost leaf in its subtree
        int range_start = left_range.first;
        int range_end = right_range.second;

        ranges[root_val] = {range_start, range_end};
        return {range_start, range_end};
    }
};

// Build the example tree from your professor
RangeTree buildExampleTree() {
    RangeTree tree;

    // Building the example tree from your professor:
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

    tree.buildWithLeaves(preorder, inorder);

    return tree;
}

int main() {
    std::cout << "Building example tree...\n\n";

    RangeTree tree = buildExampleTree();
    tree.print();

    std::cout << "\nTotal leaves: " << tree.leafCount() << std::endl;

    // Test navigation functions
    std::cout << "\nTesting navigation:\n";
    std::cout << "Parent of 1: " << tree.getParent(1) << std::endl;
    std::cout << "Left child of 4: " << tree.getLeftChild(4) << std::endl;
    std::cout << "Right child of 4: " << tree.getRightChild(4) << std::endl;

    return 0;
}