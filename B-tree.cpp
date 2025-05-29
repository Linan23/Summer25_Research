/*
Step1: Start Program w a tree notation, binary tree, how do we want to represent
Structure or array(keep track of node)

Step 2: Do any rotations on any edges

Goal: Given two sets of trees what’s the minimum of flips to get to object and how they are corresponding

When a edge is common between beginning and end you don’t flip
 - Add a leaf for each node to make a complete tree, label leaf, traverse each node
 -can directly check if they are equivalent

If we can determine a edge is same base on start and end without leafs

Write Question: Consider the binary rotation that allows me to quickly flip a edge and find similarities of edges between two trees

Brute Force as baseline, try to rotate them test from there
*/

/*
A node structure version of b tree....

1. Node def:

- each node has: val, left/righ pointers to child, pointer to parent node, indices of leaves under current node

2. Tree building logic:

two traversal arrays: preorder: root-left-right, inorder: left-root-right

It recursively:

 -Picks the current root from preorder.

 -Finds its position in inorder (using a fast hash map).

 -Splits into left and right subtrees.

 -Builds those subtrees and connects child pointers.

 -Tracks real leaves and assigns their index-based range.

3. Leaf Tracking:

If a node has no left or right child, it’s a leaf:

The code adds it to a leaves vector, then assigns its range_start and range_end to that index

Time Complexity: O(n), n is the number of nodes
Space Complexity: O(n): 
 - O(n) for nodes
 - O(n) for leaves
 - O(n) for node lookup make
 - for inorder index map
*/

#include <iostream>
#include <vector>
#include <unordered_map>

// Node definition for linked structure
struct Node {
    int val;
    Node* left;
    Node* right;
    Node* parent;

    int range_start; // Index of leftmost leaf in subtree
    int range_end;   // Index of rightmost leaf in subtree

    Node(int value) : val(value), left(nullptr), right(nullptr), parent(nullptr), range_start(-1), range_end(-1) {}
};

class RangeTreeLinked {
public:
    Node* root;
    std::vector<Node*> leaves;
    std::unordered_map<int, Node*> nodes; // quick node lookup by value

    RangeTreeLinked() : root(nullptr) {}

    // Public method to build the tree from preorder and inorder
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        leaves.clear();
        nodes.clear();
        int leaf_counter = 0;

        // Preprocess inorder into a map for O(1) lookup
        std::unordered_map<int, int> inorder_index;
        for (int i = 0; i < inorder.size(); ++i) {
            inorder_index[inorder[i]] = i;
        }

        root = buildRecursive(preorder, inorder, inorder_index,
                              0, preorder.size() - 1, 0, inorder.size() - 1, nullptr, leaf_counter);
    }

    // Print tree structure
    void printTreeStructure(Node* node, std::string indent = "") {
        if (!node) return;
        std::cout << indent << "Node " << node->val << " (range: " << node->range_start << ", " << node->range_end << ")\n";
        if (node->left || node->right) {
            printTreeStructure(node->left, indent + "  L-");
            printTreeStructure(node->right, indent + "  R-");
        }
    }

    void printLeaves() {
        std::cout << "Leaves (left to right): ";
        for (size_t i = 0; i < leaves.size(); ++i)
            std::cout << leaves[i]->val << " ";
        std::cout << std::endl;
    }

private:
    Node* buildRecursive(const std::vector<int>& preorder, const std::vector<int>& inorder,
                         const std::unordered_map<int, int>& inorder_index,
                         int pre_start, int pre_end, int in_start, int in_end,
                         Node* parent, int& leaf_counter) {

        if (pre_start > pre_end || in_start > in_end) {
            return nullptr;
        }

        int root_val = preorder[pre_start];
        Node* node = new Node(root_val);
        node->parent = parent;
        nodes[root_val] = node;

        int root_idx = inorder_index.at(root_val);
        int left_size = root_idx - in_start;

        node->left = buildRecursive(preorder, inorder, inorder_index,
                                    pre_start + 1, pre_start + left_size, in_start, root_idx - 1, node, leaf_counter);

        node->right = buildRecursive(preorder, inorder, inorder_index,
                                     pre_start + left_size + 1, pre_end, root_idx + 1, in_end, node, leaf_counter);

        if (!node->left && !node->right) {
            node->range_start = node->range_end = leaf_counter;
            leaves.push_back(node);
            leaf_counter++;
        } else {
            node->range_start = node->left ? node->left->range_start : node->right->range_start;
            node->range_end = node->right ? node->right->range_end : node->left->range_end;
        }

        return node;
    }
};

int main() {
    RangeTreeLinked tree;

    // Ge Xia example
    std::vector<int> preorder = {6, 4, 3, 2, 1, 5};
    std::vector<int> inorder = {1, 2, 3, 4, 5, 6};

    tree.build(preorder, inorder);

    std::cout << "Linked-Node Range Tree Structure:\n";
    tree.printTreeStructure(tree.root);

    tree.printLeaves();

    return 0;
}
