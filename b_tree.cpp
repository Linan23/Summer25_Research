#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

// Simple hash and equality for std::pair<int,int>
struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const {
        // Combine two 32-bit ints into one 64-bit for hashing
        return std::hash<long long>()(((long long)p.first << 32) ^ (unsigned long long)p.second);
    }
};
struct PairEq {
    bool operator()(const std::pair<int,int>& a, const std::pair<int,int>& b) const {
        return a.first == b.first && a.second == b.second;
    }
};

// Node holds a value, pointers, and the leaf-range indices.
struct Node {
    int val;
    Node* left;
    Node* right;
    Node* parent;
    int range_start;
    int range_end;

    Node(int value)
        : val(value), left(nullptr), right(nullptr), parent(nullptr),
          range_start(-1), range_end(-1) {}
};

class RangeTreeLinked {
public:
    Node* root = nullptr;
    std::vector<Node*> leaves;
    std::unordered_map<int, Node*> nodes;  // map value -> node

    // Build from preorder and inorder lists
    void build(const std::vector<int>& preorder,
               const std::vector<int>& inorder) {
        // Clear previous data (leaves and nodes)
        clearTree(root);
        leaves.clear();
        nodes.clear();
        int leaf_counter = 0;

        // map each inorder value to its index
        std::unordered_map<int,int> inorder_index;
        for(int i=0; i<(int)inorder.size(); ++i)
            inorder_index[inorder[i]] = i;

        root = buildRecursive(preorder, inorder, inorder_index,
                              0, (int)preorder.size()-1,
                              0, (int)inorder.size()-1,
                              nullptr, leaf_counter);
    }

    // Print each node and its range
    void printTreeStructure(Node* node, std::string indent = "") const {
        if(!node) return;
        std::cout << indent
                  << "Node " << node->val
                  << " (range: " << node->range_start
                  << ", " << node->range_end << ")\n";
        printTreeStructure(node->left,  indent + "  L-");
        printTreeStructure(node->right, indent + "  R-");
    }

    // Print just the leaf values in order
    void printLeaves() const {
        std::cout << "Leaves (left to right): ";
        for(auto* leaf: leaves)
            std::cout << leaf->val << " ";
        std::cout << "\n";
    }

    // Rotate x and its right child up (left-rotate)
    void rotateLeft(Node* x) {
        if(!x) return;
        Node* y = x->right;
        if(!y) return;  // no right child
        x->right = y->left;
        if(y->left) y->left->parent = x;
        y->parent = x->parent;
        if(!x->parent)                root = y;
        else if(x == x->parent->left) x->parent->left  = y;
        else                          x->parent->right = y;
        y->left = x;
        x->parent = y;
        // fix ranges up from y
        updateRangesUp(y);
    }

    // Rotate x and its left child up (right-rotate)
    void rotateRight(Node* x) {
        if(!x) return;
        Node* y = x->left;
        if(!y) return;  // no left child
        x->left = y->right;
        if(y->right) y->right->parent = x;
        y->parent = x->parent;
        if(!x->parent)                 root = y;
        else if(x == x->parent->left)  x->parent->left  = y;
        else                            x->parent->right = y;
        y->right = x;
        x->parent = y;
        // fix ranges up from y
        updateRangesUp(y);
    }

    // Collect all parent->child edges into `out_set`
    void collectEdges(Node* n,
                      std::unordered_set<std::pair<int,int>, PairHash, PairEq>& out_set) const
    {
        if(!n) return;
        if(n->left) {
            out_set.insert({n->val, n->left->val});
            collectEdges(n->left, out_set);
        }
        if(n->right) {
            out_set.insert({n->val, n->right->val});
            collectEdges(n->right, out_set);
        }
    }

    // Destructor to free all nodes
    ~RangeTreeLinked() {
        clearTree(root);
    }

private:
    // Recursively build tree and set ranges
    Node* buildRecursive(const std::vector<int>& preorder,
                         const std::vector<int>& inorder,
                         const std::unordered_map<int,int>& inorder_index,
                         int pre_start, int pre_end,
                         int in_start,  int in_end,
                         Node* parent, int& leaf_counter)
    {
        if(pre_start > pre_end || in_start > in_end)
            return nullptr;
        int root_val = preorder[pre_start];
        Node* node = new Node(root_val);
        node->parent = parent;
        nodes[root_val] = node;
        int idx = inorder_index.at(root_val);
        int left_size = idx - in_start;
        node->left = buildRecursive(preorder, inorder, inorder_index,
                                    pre_start+1, pre_start+left_size,
                                    in_start, idx-1,
                                    node, leaf_counter);
        node->right = buildRecursive(preorder, inorder, inorder_index,
                                     pre_start+left_size+1, pre_end,
                                     idx+1, in_end,
                                     node, leaf_counter);
        if(!node->left && !node->right) {
            node->range_start = node->range_end = leaf_counter;
            leaves.push_back(node);
            leaf_counter++;
        } else {
            node->range_start = node->left  ? node->left->range_start  : node->right->range_start;
            node->range_end   = node->right ? node->right->range_end : node->left->range_end;
        }
        return node;
    }

    // Update one node's range from its children
    void updateRange(Node* n) {
        if(!n->left && !n->right) return;
        n->range_start = n->left  ? n->left->range_start  : n->right->range_start;
        n->range_end   = n->right ? n->right->range_end : n->left->range_end;
    }

    // Walk up from n to root, fixing ranges
    void updateRangesUp(Node* n) {
        while(n) {
            updateRange(n);
            n = n->parent;
        }
    }

    // Recursively delete all nodes
    void clearTree(Node* n) {
        if(!n) return;
        clearTree(n->left);
        clearTree(n->right);
        delete n;
    }
};

// Compute and print shared edges between two RangeTreeLinked trees
void getSharedEdges(const RangeTreeLinked& A, const RangeTreeLinked& B) {
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
    RangeTreeLinked tree;

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
        tree.printTreeStructure(tree.root); 
        tree.printLeaves(); 
        std::cout << "\n";
        Node* x = tree.nodes[tc.rotateAt];
        if(tc.doLeft) tree.rotateLeft(x); else tree.rotateRight(x);
        std::cout << "After:\n";  
        tree.printTreeStructure(tree.root); 
        tree.printLeaves(); 
        std::cout << "\n\n";
    }

    // Test 4: Multiple rotations on 7-node
    std::cout << "Test 4: Multiple rotations on 7-node\n";
    tree.build({4,2,1,3,6,5,7}, {1,2,3,4,5,6,7});
    std::cout << "Original structure:\n"; 
    tree.printTreeStructure(tree.root); 
    tree.printLeaves(); 
    std::cout << "\n";
    tree.rotateLeft(tree.nodes[2]);
    tree.rotateRight(tree.nodes[6]);
    std::cout << "After rotateLeft(2) and rotateRight(6):\n";
    tree.printTreeStructure(tree.root); 
    tree.printLeaves(); 
    std::cout << "\n";

    // Test 5: rotateLeft on node with no right child
    std::cout << "Test 5: rotateLeft on leaf node (1)\n";
    tree.build({2,1,3}, {1,2,3});
    std::cout << "Before (should be unchanged):\n";
    tree.printTreeStructure(tree.root); 
    tree.printLeaves(); 
    std::cout << "\n";
    tree.rotateLeft(tree.nodes[1]); // 1 has no right child
    std::cout << "After (should still be unchanged):\n";
    tree.printTreeStructure(tree.root); 
    tree.printLeaves(); 
    std::cout << "\n";

    // Test 6: Shared-edge detection between two trees
    std::cout << "Test 6: Shared edges between two trees\n";
    // Build a second tree to compare
    RangeTreeLinked tree2;
    tree.build({4,2,6,5,7}, {2,4,5,6,7}); // Tree1
    tree2.build({4,2,6,3,5}, {2,4,3,6,5}); // Tree2
    std::cout << "Tree1 structure:\n";
    tree.printTreeStructure(tree.root);
    tree.printLeaves();
    std::cout << "\nTree2 structure:\n";
    tree2.printTreeStructure(tree2.root);
    tree2.printLeaves();
    std::cout << "\nShared edges:\n";
    getSharedEdges(tree, tree2);

    std::cout << "\n Additional Tests for shared edge\n\n";

    // Both trees empty
    {
        std::cout << "-- Case 1: Both trees empty --\n";
        RangeTreeLinked A, B;
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
        RangeTreeLinked A, B;
        A.build({}, {});                    // A is empty
        B.build({1,2,3}, {2,1,3});          // B:   1
                                            //      / \
                                            //     2   3
        std::cout << "Tree A (empty), Tree B:\n";
        std::cout << "Tree B structure:\n";
        B.printTreeStructure(B.root);
        B.printLeaves();
        std::cout << "Shared edges (expected none):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // Identical small tree
    {
        std::cout << "Case 3: Identical small tree\n";
        RangeTreeLinked A, B;
        A.build({2,1,3}, {1,2,3});  //      2
                                    //     / \
                                    //    1   3
        B.build({2,1,3}, {1,2,3});
        std::cout << "Tree A and Tree B (identical):\n";
        A.printTreeStructure(A.root);
        A.printLeaves();
        std::cout << "Shared edges (expected: 2->1, 2->3):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // Same values, different shape
    {
        std::cout << "Case 4: Same values, different shape\n";
        RangeTreeLinked A, B;
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
        A.printTreeStructure(A.root);
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
        B.printTreeStructure(B.root);
        B.printLeaves();

        std::cout << "Shared edges (expected none):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }

    // Partial overlap deep inside
    {
        std::cout << "Case 5: Partial overlap deep inside\n";
        RangeTreeLinked A, B;
        A.build({5,3,1,8,7,9}, {1,3,5,7,8,9});
        // Tree A:
        //      5
        //     / \
        //    3   8
        //   /   / \
        //  1   7   9
        std::cout << "Tree A structure:\n";
        A.printTreeStructure(A.root);
        A.printLeaves();

        B.build({5,2,8,7}, {2,5,7,8});
        // Tree B:
        //      5
        //     / \
        //    2   8
        //       /
        //      7
        std::cout << "Tree B structure:\n";
        B.printTreeStructure(B.root);
        B.printLeaves();

        std::cout << "Shared edges (expected: 5->8, 8->7):\n";
        getSharedEdges(A, B);
        std::cout << "\n";
    }


    return 0;
}
