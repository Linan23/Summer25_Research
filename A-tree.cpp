/**
 * How A-tree works:
 *
 * Step 1: Leaf placement
 *      - Traverse the tree recursively, whenever a node is missing a child -> add a leaf node
 * then number leaves from left to right as they are encountered.
 *
 * Step 2: Build Leaf Array
 *      - leaves vector stores all leaf values in left-to-right order, each leaf get a pos index
 *
 * Step 3: Calculate Ranges:
 *      - Each node range is the (leftmost left pos, rightmost left pos) in its subtree
 *      - Range represents leaf span
 */
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>

struct RangeTree {
    std::vector<int> leaves;                   // Store leaf values in left-to-right order
    std::unordered_map<int, std::pair<int, int>> ranges; // node_value -> (start, end) leaf range

    RangeTree() {}

    // Build tree and calculate ranges based on leaf positions
    void buildWithLeaves(const std::vector<int>& preorder, const std::vector<int>& inorder) {
        leaves.clear();
        ranges.clear();
        if (!preorder.empty()) {
            buildRecursiveWithLeaves(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1);
        }
    }

    // Print the tree representation
    void print() const {
        std::cout << "Leaves in left-to-right order: ";
        for (size_t i = 0; i < leaves.size(); i++) {
            std::cout << "L" << i << " ";
        }
        std::cout << "\nRanges (actual nodes only):\n";

        // Print in the same order as the original tree structure
        std::vector<int> node_order = {1, 2, 3, 4, 5, 6}; // Adjust based on your tree
        for (int node : node_order) {
            if (ranges.find(node) != ranges.end()) {
                auto range = ranges.at(node);
                std::cout << node << "(" << range.first << "," << range.second << ") ";
            }
        }
        std::cout << std::endl;
    }

    // Get the total number of leaves
    int leafCount() const {
        return leaves.size();
    }

private:
    std::pair<int, int> buildRecursiveWithLeaves(const std::vector<int>& preorder, const std::vector<int>& inorder,
                                                 int pre_start, int pre_end, int in_start, int in_end) {
        if (pre_start > pre_end || in_start > in_end) {
            // Add a leaf and return its position
            int leaf_pos = leaves.size();
            leaves.push_back(1000 + leaf_pos); // Unique leaf value
            return {leaf_pos, leaf_pos};
        }

        int root_val = preorder[pre_start];

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
        std::pair<int, int> left_range = buildRecursiveWithLeaves(preorder, inorder,
                                                                  pre_start + 1, pre_start + left_size,
                                                                  in_start, root_idx - 1);

        // Build right subtree
        std::pair<int, int> right_range = buildRecursiveWithLeaves(preorder, inorder,
                                                                   pre_start + left_size + 1, pre_end,
                                                                   root_idx + 1, in_end);

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
    //        /
    //       4
    //      /\
    //     3  5
    //    /
    //   2
    //  /
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

    return 0;
}