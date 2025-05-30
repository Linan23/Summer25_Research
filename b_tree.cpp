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

// Node holds a value, pointers, and the leaf-range indices.
struct Node
{
    int val;
    Node *left;
    Node *right;
    Node *parent;
    int range_start;
    int range_end;

    Node(int value)
        : val(value), left(nullptr), right(nullptr), parent(nullptr),
          range_start(-1), range_end(-1) {}
};

class RangeTreeLinked
{
public:
    Node *root = nullptr;
    std::vector<Node *> leaves;
    std::unordered_map<int, Node *> nodes; // map value -> node

    // Build from preorder and inorder lists
    void build(const std::vector<int> &preorder,
               const std::vector<int> &inorder)
    {
        leaves.clear();
        nodes.clear();
        int leaf_counter = 0;

        // map each inorder value to its index
        std::unordered_map<int, int> inorder_index;
        for (int i = 0; i < (int)inorder.size(); ++i)
        {
            inorder_index[inorder[i]] = i;
        }

        root = buildRecursive(preorder, inorder, inorder_index,
                              0, (int)preorder.size() - 1,
                              0, (int)inorder.size() - 1,
                              nullptr, leaf_counter);
    }

    // Print each node and its range
    void printTreeStructure(Node *node, std::string indent = "")
    {
        if (!node)
            return;
        std::cout << indent
                  << "Node " << node->val
                  << " (range: " << node->range_start
                  << ", " << node->range_end << ")\n";
        printTreeStructure(node->left, indent + "  L-");
        printTreeStructure(node->right, indent + "  R-");
    }

    // Print just the leaf values in order
    void printLeaves()
    {
        std::cout << "Leaves (left to right): ";
        for (auto *leaf : leaves)
            std::cout << leaf->val << " ";
        std::cout << "\n";
    }

    // Rotate x and its right child up (left-rotate)
    void rotateLeft(Node *x)
    {
        Node *y = x->right;
        if (!y)
            return; // no right child

        // move y's left subtree to x.right
        x->right = y->left;
        if (y->left)
            y->left->parent = x;

        // link y into x's old parent
        y->parent = x->parent;
        if (!x->parent)
            root = y;
        else if (x == x->parent->left)
            x->parent->left = y;
        else
            x->parent->right = y;

        // finish rotation
        y->left = x;
        x->parent = y;

        // fix ranges up from y
        updateRangesUp(y);
    }

    // Rotate x and its left child up (right-rotate)
    void rotateRight(Node *x)
    {
        Node *y = x->left;
        if (!y)
            return; // no left child

        // move y's right subtree to x.left
        x->left = y->right;
        if (y->right)
            y->right->parent = x;

        // link y into x's old parent
        y->parent = x->parent;
        if (!x->parent)
            root = y;
        else if (x == x->parent->left)
            x->parent->left = y;
        else
            x->parent->right = y;

        // finish rotation
        y->right = x;
        x->parent = y;

        // fix ranges up from y
        updateRangesUp(y);
    }

private:
    // Recursively build tree and set ranges
    Node *buildRecursive(const std::vector<int> &preorder,
                         const std::vector<int> &inorder,
                         const std::unordered_map<int, int> &inorder_index,
                         int pre_start, int pre_end,
                         int in_start, int in_end,
                         Node *parent, int &leaf_counter)
    {
        if (pre_start > pre_end || in_start > in_end)
            return nullptr;

        int root_val = preorder[pre_start];
        Node *node = new Node(root_val);
        node->parent = parent;
        nodes[root_val] = node;

        int idx = inorder_index.at(root_val);
        int left_size = idx - in_start;

        // build left subtree
        node->left = buildRecursive(preorder, inorder, inorder_index,
                                    pre_start + 1,
                                    pre_start + left_size,
                                    in_start,
                                    idx - 1,
                                    node,
                                    leaf_counter);

        // build right subtree
        node->right = buildRecursive(preorder, inorder, inorder_index,
                                     pre_start + left_size + 1,
                                     pre_end,
                                     idx + 1,
                                     in_end,
                                     node,
                                     leaf_counter);

        // set range
        if (!node->left && !node->right)
        {
            // leaf
            node->range_start = node->range_end = leaf_counter;
            leaves.push_back(node);
            leaf_counter++;
        }
        else
        {
            node->range_start = node->left ? node->left->range_start
                                           : node->right->range_start;
            node->range_end = node->right ? node->right->range_end
                                          : node->left->range_end;
        }

        return node;
    }

    // Update one node's range from its children
    void updateRange(Node *n)
    {
        if (!n->left && !n->right)
            return; // leaf keeps its own
        n->range_start = n->left ? n->left->range_start
                                 : n->right->range_start;
        n->range_end = n->right ? n->right->range_end
                                : n->left->range_end;
    }

    // Walk up from n to root, fixing ranges
    void updateRangesUp(Node *n)
    {
        while (n)
        {
            updateRange(n);
            n = n->parent;
        }
    }
};

int main()
{
    RangeTreeLinked tree;

    struct TestCase
    {
        std::string name;
        std::vector<int> preorder;
        std::vector<int> inorder;
        int rotateAt;
        bool doLeft; // true = rotateLeft, false = rotateRight
    };

    std::vector<TestCase> tests = {
        {"Test 1: Simple 3-node, left-rotate at 2",
         {2, 1, 3},
         {1, 2, 3},
         2,
         true},
        {"Test 2: Right-chain of 4, left-rotate at 2",
         {1, 2, 3, 4},
         {1, 2, 3, 4},
         2,
         true},
        {"Test 3: Balanced 7-node, right-rotate at 6",
         {4, 2, 1, 3, 6, 5, 7},
         {1, 2, 3, 4, 5, 6, 7},
         6,
         false}

    };

    for (auto &tc : tests)
    {
        std::cout << "=== " << tc.name << " ===\n";
        tree.build(tc.preorder, tc.inorder);

        std::cout << "Before:\n";
        tree.printTreeStructure(tree.root);
        tree.printLeaves();
        std::cout << "\n";

        Node *x = tree.nodes[tc.rotateAt];
        if (tc.doLeft)
            tree.rotateLeft(x);
        else
            tree.rotateRight(x);

        std::cout << "After:\n";
        tree.printTreeStructure(tree.root);
        tree.printLeaves();
        std::cout << "\n\n";
    }

    // --- Test 4: Multiple rotations on 7-node ---
    std::cout << "=== Test 4: Multiple rotations on 7-node ===\n";
    tree.build(
        /* preorder */ {4, 2, 1, 3, 6, 5, 7},
        /* inorder  */ {1, 2, 3, 4, 5, 6, 7});

    std::cout << "Original structure:\n";
    tree.printTreeStructure(tree.root);
    tree.printLeaves();
    std::cout << "\n";

    // First, left-rotate at node 2:
    tree.rotateLeft(tree.nodes[2]);
    // Then, right-rotate at node 6:
    tree.rotateRight(tree.nodes[6]);

    std::cout << "After rotateLeft(2) and rotateRight(6):\n";
    tree.printTreeStructure(tree.root);
    tree.printLeaves();
    std::cout << "\n";
    
    return 0;
}