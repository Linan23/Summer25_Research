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
    // Vector-based storage for MASSIVE performance improvement - O(1) access instead of O(log n) hash lookups
    std::vector<std::pair<int, int>> ranges;  // index = node_id, value = (start, end) range
    std::vector<int> parent;                  // index = node_id, value = parent_id
    std::vector<int> left_child;              // index = node_id, value = left_child_id
    std::vector<int> right_child;             // index = node_id, value = right_child_id
    std::vector<bool> is_dummy;               // index = node_id, value = true if dummy leaf
    std::vector<bool> node_exists;            // index = node_id, value = true if node exists
    std::vector<int> position_in_inorder;     // index = node_id, value = position in original inorder

    int root;                                 // root node ID
    std::vector<int> leaves;                  // leaf node IDs in left-to-right order

    // Store the original inorder sequence for range calculations
    std::vector<int> original_inorder;

    // Mapping from original values to internal IDs (needed since original values might not be 0-based)
    std::unordered_map<int, int> value_to_id; // original value -> internal ID
    std::unordered_map<int, int> id_to_value; // internal ID -> original value
    int next_id;                              // next available ID for original nodes
    int next_dummy_id;                        // next available ID for dummy nodes (negative)

    // Constants for "no value" - cleaner than magic numbers
    static constexpr int NO_PARENT = -1;
    static constexpr int NO_CHILD = -1;

    VectorRangeTreeMap() : root(NO_PARENT), next_id(0), next_dummy_id(-1000) {}

    // Build from preorder and inorder lists
    void build(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        // Clear previous data (leaves and nodes)
        clear();

        if (preorder.empty()) return;

        // Create ID mapping for arbitrary values - converts user values to 0-based IDs
        createIdMapping(preorder, inorder);

        // Convert to ID-based sequences for internal processing
        std::vector<int> preorder_ids, inorder_ids;
        convertToIds(preorder, preorder_ids);
        convertToIds(inorder, inorder_ids);

        // Store original inorder as IDs
        original_inorder = inorder_ids;

        // Create position mapping from ORIGINAL inorder sequence - O(n) preprocessing for O(1) lookups
        for (size_t i = 0; i < inorder_ids.size(); i++) {
            position_in_inorder[inorder_ids[i]] = i;
        }

        root = preorder_ids[0];
        buildRecursive(preorder_ids, inorder_ids, 0, preorder_ids.size() - 1, 0, inorder_ids.size() - 1, NO_PARENT);

        // Build leaves vector in left-to-right order
        buildLeavesVector(root);

        addDummyLeaves();
        calculateAllRanges();
    }

    // ULTRA-FAST lookups - single array access, no hashing! O(1) guaranteed
    inline int getParent(int node_id) const {
        int index = (node_id >= 0) ? node_id : -node_id;
        return (index < parent.size()) ? parent[index] : NO_PARENT;
    }

    inline int getLeftChild(int node_id) const {
        int index = (node_id >= 0) ? node_id : -node_id;
        return (index < left_child.size()) ? left_child[index] : NO_CHILD;
    }

    inline int getRightChild(int node_id) const {
        int index = (node_id >= 0) ? node_id : -node_id;
        return (index < right_child.size()) ? right_child[index] : NO_CHILD;
    }

    inline bool isDummy(int node_id) const {
        int index = (node_id >= 0) ? node_id : -node_id;
        return (index < is_dummy.size()) ? is_dummy[index] : false;
    }

    inline bool exists(int node_id) const {
        int index = (node_id >= 0) ? node_id : -node_id;
        return (index < node_exists.size()) ? node_exists[index] : false;
    }

    inline std::pair<int, int> getRange(int node_id) const {
        int index = (node_id >= 0) ? node_id : -node_id;
        return (index < ranges.size()) ? ranges[index] : std::make_pair(0, 0);
    }

    // Print each node and its range
    void printTreeStructure(int node_id = -2, std::string indent = "") const {
        if (node_id == -2) node_id = root; // Default to root on first call
        if (node_id == NO_PARENT || !exists(node_id)) return;

        // Convert back to original value for display
        int original_value = id_to_value.at(node_id);
        auto range = getRange(node_id);

        std::cout << indent << "Node " << original_value
                  << " (range: " << range.first << ", " << range.second << ")\n";

        int left = getLeftChild(node_id);
        if (left != NO_CHILD && !isDummy(left)) {
            printTreeStructure(left, indent + "  L-");
        }

        int right = getRightChild(node_id);
        if (right != NO_CHILD && !isDummy(right)) {
            printTreeStructure(right, indent + "  R-");
        }
    }

    // Print just the leaf values in order
    void printLeaves() const {
        std::cout << "Leaves (left to right): ";
        for (int leaf_id : leaves) {
            if (id_to_value.find(leaf_id) != id_to_value.end()) {
                std::cout << id_to_value.at(leaf_id) << " ";
            }
        }
        std::cout << "\n";
    }

    // Rotate x and its right child up (left-rotate) - Only updates the two affected nodes
    void rotateLeft(int x_value) {
        // Convert from original value to internal ID
        if (value_to_id.find(x_value) == value_to_id.end()) return;
        int x = value_to_id[x_value];

        if (x == NO_PARENT || !exists(x)) return;

        int y = getRightChild(x);  // Single array access!
        if (y == NO_CHILD || isDummy(y)) return; // can't rotate dummy

        // Perform rotation logic
        int y_left = getLeftChild(y);
        int x_parent = getParent(x);

        // Update x's right child
        int x_index = (x >= 0) ? x : -x;
        if (y_left != NO_CHILD) {
            right_child[x_index] = y_left;
            int y_left_index = (y_left >= 0) ? y_left : -y_left;
            parent[y_left_index] = x;
        } else {
            right_child[x_index] = NO_CHILD;
        }

        // Update y's parent
        int y_index = (y >= 0) ? y : -y;
        if (x_parent != NO_PARENT) {
            parent[y_index] = x_parent;
            int x_parent_index = (x_parent >= 0) ? x_parent : -x_parent;

            if (getLeftChild(x_parent) == x) {
                left_child[x_parent_index] = y;
            } else {
                right_child[x_parent_index] = y;
            }
        } else {
            root = y;
            parent[y_index] = NO_PARENT;
        }

        // Set y as parent of x
        left_child[y_index] = x;
        parent[x_index] = y;

        // FIXED: Only update ranges for the two nodes involved in rotation
        updateNodeRange(x);  // Update x's range first (now child)
        updateNodeRange(y);  // Update y's range second (now parent)
    }

    // Rotate x and its left child up (right-rotate) - Only updates the two affected nodes
    void rotateRight(int x_value) {
        // Convert from original value to internal ID
        if (value_to_id.find(x_value) == value_to_id.end()) return;
        int x = value_to_id[x_value];

        if (x == NO_PARENT || !exists(x)) return;

        int y = getLeftChild(x);
        if (y == NO_CHILD || isDummy(y)) return; // can't rotate dummy

        // Perform rotation logic
        int y_right = getRightChild(y);
        int x_parent = getParent(x);

        // Update x's left child
        int x_index = (x >= 0) ? x : -x;
        if (y_right != NO_CHILD) {
            left_child[x_index] = y_right;
            int y_right_index = (y_right >= 0) ? y_right : -y_right;
            parent[y_right_index] = x;
        } else {
            left_child[x_index] = NO_CHILD;
        }

        // Update y's parent
        int y_index = (y >= 0) ? y : -y;
        if (x_parent != NO_PARENT) {
            parent[y_index] = x_parent;
            int x_parent_index = (x_parent >= 0) ? x_parent : -x_parent;

            if (getLeftChild(x_parent) == x) {
                left_child[x_parent_index] = y;
            } else {
                right_child[x_parent_index] = y;
            }
        } else {
            root = y;
            parent[y_index] = NO_PARENT;
        }

        // Set y as parent of x
        right_child[y_index] = x;
        parent[x_index] = y;

        // FIXED: Only update ranges for the two nodes involved in rotation
        updateNodeRange(x);  // Update x's range first (now child)
        updateNodeRange(y);  // Update y's range second (now parent)
    }

    // Collect all parent->child edges into `out_set`
    void collectEdges(int node_id,
                      std::unordered_set<std::pair<int,int>, PairHash, PairEq>& out_set) const
    {
        if (node_id == NO_PARENT || !exists(node_id)) return;

        // Convert back to original values for edge collection
        int original_node = id_to_value.at(node_id);

        int left = getLeftChild(node_id);
        if (left != NO_CHILD && !isDummy(left)) {
            int original_left = id_to_value.at(left);
            out_set.insert({original_node, original_left});
            collectEdges(left, out_set);
        }

        int right = getRightChild(node_id);
        if (right != NO_CHILD && !isDummy(right)) {
            int original_right = id_to_value.at(right);
            out_set.insert({original_node, original_right});
            collectEdges(right, out_set);
        }
    }

    // Print ranges for all original nodes
    void print() const {
        std::cout << "Tree structure:\n";
        printTreeStructure();

        std::cout << "\nRanges (original nodes only):\n";

        // Get all original nodes and sort them
        std::vector<int> original_values;
        for (const auto& [id, value] : id_to_value) {
            if (!isDummy(id)) {
                original_values.push_back(value);
            }
        }
        std::sort(original_values.begin(), original_values.end());

        // Print all original nodes
        for (int value : original_values) {
            int id = value_to_id.at(value);
            auto range = getRange(id);
            std::cout << value << "(" << range.first << "," << range.second << ") ";
        }
        std::cout << std::endl;
    }

private:
    void clear() {
        ranges.clear();
        parent.clear();
        left_child.clear();
        right_child.clear();
        is_dummy.clear();
        node_exists.clear();
        position_in_inorder.clear();
        leaves.clear();
        original_inorder.clear();
        value_to_id.clear();
        id_to_value.clear();
        root = NO_PARENT;
        next_id = 0;
        next_dummy_id = -1000;
    }

    // Create bidirectional mapping between original values and internal IDs
    void createIdMapping(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        std::set<int> unique_values;
        for (int val : preorder) unique_values.insert(val);
        for (int val : inorder) unique_values.insert(val);

        for (int value : unique_values) {
            int id = next_id++;
            value_to_id[value] = id;
            id_to_value[id] = value;
        }

        // Resize vectors to accommodate all IDs plus dummy nodes
        int max_size = next_id + 1000; // Extra space for dummy nodes
        ranges.resize(max_size, {0, 0});
        parent.resize(max_size, NO_PARENT);
        left_child.resize(max_size, NO_CHILD);
        right_child.resize(max_size, NO_CHILD);
        is_dummy.resize(max_size, false);
        node_exists.resize(max_size, false);
        position_in_inorder.resize(max_size, -1);

        // Mark original nodes as existing
        for (const auto& [value, id] : value_to_id) {
            node_exists[id] = true;
        }
    }

    void convertToIds(const std::vector<int>& values, std::vector<int>& ids) {
        ids.clear();
        ids.reserve(values.size());
        for (int value : values) {
            ids.push_back(value_to_id[value]);
        }
    }

    // Build tree recursively
    void buildRecursive(const std::vector<int>& preorder_ids, const std::vector<int>& inorder_ids,
                        int pre_start, int pre_end, int in_start, int in_end, int parent_id) {
        if (pre_start > pre_end || in_start > in_end) return;

        int root_id = preorder_ids[pre_start];

        // Set parent relationship
        if (parent_id != NO_PARENT) {
            int root_index = (root_id >= 0) ? root_id : -root_id;
            parent[root_index] = parent_id;
        }

        // Find root position in inorder - O(1) lookup instead of O(n) search!
        int root_idx = position_in_inorder[root_id];
        int left_size = root_idx - in_start;

        // Build left subtree
        if (left_size > 0) {
            int left_root_id = preorder_ids[pre_start + 1];
            int root_index = (root_id >= 0) ? root_id : -root_id;
            left_child[root_index] = left_root_id;
            buildRecursive(preorder_ids, inorder_ids, pre_start + 1, pre_start + left_size,
                           in_start, root_idx - 1, root_id);
        }

        // Build right subtree
        if (root_idx < in_end) {
            int right_root_id = preorder_ids[pre_start + left_size + 1];
            int root_index = (root_id >= 0) ? root_id : -root_id;
            right_child[root_index] = right_root_id;
            buildRecursive(preorder_ids, inorder_ids, pre_start + left_size + 1, pre_end,
                           root_idx + 1, in_end, root_id);
        }
    }

    // Build leaves vector in left-to-right order
    void buildLeavesVector(int node_id) {
        if (node_id == NO_PARENT || !exists(node_id)) return;

        int left = getLeftChild(node_id);
        int right = getRightChild(node_id);

        // If it's a leaf (no children), add to leaves
        if (left == NO_CHILD && right == NO_CHILD) {
            leaves.push_back(node_id);
            return;
        }

        // Recursively build leaves in inorder
        if (left != NO_CHILD) buildLeavesVector(left);
        if (right != NO_CHILD) buildLeavesVector(right);
    }

    // Add dummy leaves to make every original node an inner node
    void addDummyLeaves() {
        std::vector<int> original_nodes = original_inorder;

        // Add dummy leaves to nodes that don't have both children
        for (int node_id : original_nodes) {
            int left = getLeftChild(node_id);
            int right = getRightChild(node_id);

            // Add dummy left child if missing
            if (left == NO_CHILD) {
                int dummy_id = addDummyNode();
                int node_index = (node_id >= 0) ? node_id : -node_id;
                left_child[node_index] = dummy_id;
                int dummy_index = (dummy_id >= 0) ? dummy_id : -dummy_id;
                parent[dummy_index] = node_id;
            }

            // Add dummy right child if missing
            if (right == NO_CHILD) {
                int dummy_id = addDummyNode();
                int node_index = (node_id >= 0) ? node_id : -node_id;
                right_child[node_index] = dummy_id;
                int dummy_index = (dummy_id >= 0) ? dummy_id : -dummy_id;
                parent[dummy_index] = node_id;
            }
        }

        // Recursively add dummy leaves to newly created inner nodes if needed
        addDummyLeavesRecursive();
    }

    int addDummyNode() {
        int dummy_id = next_dummy_id--;

        // Ensure vectors are large enough
        int index = -dummy_id;
        if (index >= ranges.size()) {
            int new_size = index + 100;
            ranges.resize(new_size, {0, 0});
            parent.resize(new_size, NO_PARENT);
            left_child.resize(new_size, NO_CHILD);
            right_child.resize(new_size, NO_CHILD);
            is_dummy.resize(new_size, false);
            node_exists.resize(new_size, false);
            position_in_inorder.resize(new_size, -1);
        }

        is_dummy[index] = true;
        node_exists[index] = true;

        return dummy_id;
    }

    // Recursively ensure all non-dummy nodes have two children
    void addDummyLeavesRecursive() {
        bool added_new_dummies = true;

        while (added_new_dummies) {
            added_new_dummies = false;
            std::vector<int> nodes_to_check;

            // Collect all non-dummy nodes
            for (const auto& [value, id] : value_to_id) {
                if (!isDummy(id)) {
                    nodes_to_check.push_back(id);
                }
            }

            // Add root if it exists
            if (root != NO_PARENT && !isDummy(root)) {
                nodes_to_check.push_back(root);
            }

            // Remove duplicates
            std::sort(nodes_to_check.begin(), nodes_to_check.end());
            nodes_to_check.erase(std::unique(nodes_to_check.begin(), nodes_to_check.end()), nodes_to_check.end());

            for (int node_id : nodes_to_check) {
                if (isDummy(node_id)) continue;

                int left = getLeftChild(node_id);
                int right = getRightChild(node_id);

                // Add dummy left child if missing
                if (left == NO_CHILD) {
                    int dummy_id = addDummyNode();
                    int node_index = (node_id >= 0) ? node_id : -node_id;
                    left_child[node_index] = dummy_id;
                    int dummy_index = (dummy_id >= 0) ? dummy_id : -dummy_id;
                    parent[dummy_index] = node_id;
                    added_new_dummies = true;
                }

                // Add dummy right child if missing
                if (right == NO_CHILD) {
                    int dummy_id = addDummyNode();
                    int node_index = (node_id >= 0) ? node_id : -node_id;
                    right_child[node_index] = dummy_id;
                    int dummy_index = (dummy_id >= 0) ? dummy_id : -dummy_id;
                    parent[dummy_index] = node_id;
                    added_new_dummies = true;
                }
            }
        }
    }

    // TRUE O(1) range update for a single node (with dummy leaves)
    void updateNodeRange(int node_id) {
        int index = (node_id >= 0) ? node_id : -node_id;

        if (isDummy(node_id)) {
            // Dummy leaves have empty ranges
            ranges[index] = {0, 0};
            return;
        }

        // For inner nodes with dummy leaves, range calculation is O(1)
        int left = getLeftChild(node_id);
        int right = getRightChild(node_id);

        // With dummy leaves, every original node is guaranteed to have both children
        if (left == NO_CHILD || right == NO_CHILD) {
            // This shouldn't happen with proper dummy leaf augmentation
            int pos = position_in_inorder[node_id];
            ranges[index] = {pos, pos + 1};
            return;
        }

        // O(1) calculation: combine ranges from left and right children
        auto left_range = getRange(left);
        auto right_range = getRange(right);
        int pos = position_in_inorder[node_id];

        int start = std::min(left_range.first, pos);
        int end = std::max(right_range.second, pos + 1);

        // If left child is dummy, start from this node's position
        if (isDummy(left)) {
            start = pos;
        }

        // If right child is dummy, end at this node's position + 1
        if (isDummy(right)) {
            end = pos + 1;
        }

        ranges[index] = {start, end};
    }


    // Calculate ranges for all nodes - only called during initial build
    void calculateAllRanges() {
        ranges.clear();
        ranges.resize(ranges.capacity(), {0, 0});

        // Initialize dummy leaf ranges (empty ranges)
        for (const auto& [value, id] : value_to_id) {
            if (isDummy(id)) {
                int index = (id >= 0) ? id : -id;
                ranges[index] = {0, 0};
            }
        }

        // Calculate ranges for original nodes in post-order (bottom-up)
        calculateRangesPostOrder(root);
    }

    // Post-order traversal to calculate ranges bottom-up
    void calculateRangesPostOrder(int node_id) {
        if (node_id == NO_PARENT || !exists(node_id)) return;

        int left = getLeftChild(node_id);
        int right = getRightChild(node_id);

        // Process children first (post-order)
        if (left != NO_CHILD) {
            calculateRangesPostOrder(left);
        }
        if (right != NO_CHILD) {
            calculateRangesPostOrder(right);
        }

        // Now calculate this node's range
        updateNodeRange(node_id);
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

    struct TestCase {
        std::string name;
        std::vector<int> preorder;
        std::vector<int> inorder;
        int rotateAt;
        bool doLeft;
    };

    std::vector<TestCase> tests = {
            {"Test 1: Simple 3-node, left-rotate at 2",    {2, 1, 3},             {1, 2, 3},             2, true},
            {"Test 2: Right-chain of 4, left-rotate at 2", {1, 2, 3, 4},          {1, 2, 3, 4},          2, true},
            {"Test 3: Balanced 7-node, right-rotate at 6", {4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7}, 6, false}
    };

    for (auto &tc: tests) {
        std::cout << tc.name << "\n";
        tree.build(tc.preorder, tc.inorder);
        std::cout << "Before:\n";
        tree.printTreeStructure();
        tree.printLeaves();
        std::cout << "\n";
        int x = tc.rotateAt; // Use the original value directly
        if (tc.doLeft) tree.rotateLeft(x); else tree.rotateRight(x);
        std::cout << "After:\n";
        tree.printTreeStructure();
        tree.printLeaves();
        std::cout << "\n\n";
    }

    // Test 4: Multiple rotations on 7-node
    std::cout << "Test 4: Multiple rotations on 7-node\n";
    tree.build({4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7});
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
    tree.build({2, 1, 3}, {1, 2, 3});
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
    VectorRangeTreeMap tree2;
    tree.build({4, 2, 6, 5, 7}, {2, 4, 5, 6, 7}); // Tree1
    tree2.build({4, 2, 6, 3, 5}, {2, 4, 3, 6, 5}); // Tree2
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
        VectorRangeTreeMap A, B;
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
        VectorRangeTreeMap A, B;
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
        VectorRangeTreeMap A, B;
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
        VectorRangeTreeMap A, B;
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
        VectorRangeTreeMap A, B;
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