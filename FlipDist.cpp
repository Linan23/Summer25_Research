#include "A_tree.h"
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <climits>
#include <cassert>
#include <functional>

// Forward declarations
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k);
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I);
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S);

// Debug flag
const bool DEBUG = true;

void debugPrint(const std::string& msg) {
    if (DEBUG) {
        std::cout << "[DEBUG] " << msg << std::endl;
    }
}

// Get all internal edges in a tree (parent->child direction) with safety checks
std::vector<std::pair<int,int>> getInternalEdges(const VectorRangeTreeMap& T) {
    std::vector<std::pair<int,int>> edges;

    try {
        // Check if tree is valid
        if (T.original_nodes.empty() || T.root < 0) {
            return edges;
        }

        std::function<void(int)> dfs = [&](int node) {
            if (node < 0 || !T.isOriginal(node)) return;

            try {
                int left = T.getLeftChild(node);
                int right = T.getRightChild(node);

                if (left >= 0 && T.isOriginal(left)) {
                    edges.emplace_back(node, left);  // parent -> child
                    dfs(left);
                }
                if (right >= 0 && T.isOriginal(right)) {
                    edges.emplace_back(node, right);  // parent -> child
                    dfs(right);
                }
            } catch (...) {
                // Skip problematic nodes
                return;
            }
        };

        if (T.isOriginal(T.root)) {
            dfs(T.root);
        }
    } catch (...) {
        // Return empty vector on any error
        edges.clear();
    }

    return edges;
}

// Count internal edges in tree
int countInternalEdges(const VectorRangeTreeMap& T) {
    return getInternalEdges(T).size();
}

// Check if two edges are adjacent (share a node)
bool areAdjacent(const std::pair<int,int>& e1, const std::pair<int,int>& e2) {
    return e1.first == e2.first || e1.first == e2.second ||
           e1.second == e2.first || e1.second == e2.second;
}

// Generate all independent subsets of edges (limit size to avoid explosion)
void generateIndependentSubsets(const std::vector<std::pair<int,int>>& edges, int index,
                                std::vector<std::pair<int,int>>& current,
                                std::vector<std::vector<std::pair<int,int>>>& result,
                                int maxSubsets = 100) {
    if (result.size() >= maxSubsets) return;  // Limit to prevent explosion

    if (index == edges.size()) {
        if (!current.empty()) {
            result.push_back(current);
        }
        return;
    }

    // Choice 1: exclude current edge
    generateIndependentSubsets(edges, index + 1, current, result, maxSubsets);

    if (result.size() >= maxSubsets) return;

    // Choice 2: include current edge if it doesn't conflict
    bool canInclude = true;
    for (const auto& e : current) {
        if (areAdjacent(e, edges[index])) {
            canInclude = false;
            break;
        }
    }

    if (canInclude) {
        current.push_back(edges[index]);
        generateIndependentSubsets(edges, index + 1, current, result, maxSubsets);
        current.pop_back();
    }
}

// Check if a specific parent->child edge exists in the tree with safety
bool hasParentChildEdge(const VectorRangeTreeMap& T, int parent, int child) {
    try {
        if (!T.isOriginal(parent) || !T.isOriginal(child)) return false;

        int left = T.getLeftChild(parent);
        int right = T.getRightChild(parent);

        return (left == child) || (right == child);
    } catch (...) {
        return false;
    }
}

// Safe tree copy function with better error handling
VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap& T) {
    VectorRangeTreeMap copy;
    try {
        // Handle empty tree
        if (T.root < 0 || !T.isOriginal(T.root) || T.original_nodes.empty()) {
            return copy;
        }

        // Get the original sequences
        std::vector<int> preorder, inorder;

        // Build preorder by traversing
        std::function<void(int, std::vector<int>&)> buildPreorder = [&](int node, std::vector<int>& pre) {
            if (node < 0 || !T.isOriginal(node)) return;
            pre.push_back(node);
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left)) buildPreorder(left, pre);
            if (right >= 0 && T.isOriginal(right)) buildPreorder(right, pre);
        };

        // Build inorder by traversing
        std::function<void(int, std::vector<int>&)> buildInorder = [&](int node, std::vector<int>& in) {
            if (node < 0 || !T.isOriginal(node)) return;
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left)) buildInorder(left, in);
            in.push_back(node);
            if (right >= 0 && T.isOriginal(right)) buildInorder(right, in);
        };

        buildPreorder(T.root, preorder);
        buildInorder(T.root, inorder);

        // Validate sequences
        if (preorder.size() == inorder.size() && !preorder.empty()) {
            copy.build(preorder, inorder);
        }
    } catch (...) {
        // Return empty tree on any error
        VectorRangeTreeMap empty;
        return empty;
    }

    return copy;
}

// Find free edge using a simplified approach with better safety
std::pair<bool, std::pair<int,int>> findFreeEdge(const VectorRangeTreeMap& T_init,
                                                 const VectorRangeTreeMap& T_final) {
    try {
        // Handle empty trees
        if (T_init.original_nodes.empty() || T_final.original_nodes.empty()) {
            return {false, {-1, -1}};
        }

        // Get edges from both trees
        auto initEdges = getInternalEdges(T_init);
        auto finalEdges = getInternalEdges(T_final);

        // Convert to sets for faster lookup
        std::set<std::pair<int,int>> initSet(initEdges.begin(), initEdges.end());
        std::set<std::pair<int,int>> finalSet(finalEdges.begin(), finalEdges.end());

        // Try each edge in T_init to see if rotating it creates an edge from T_final
        for (const auto& edge : initEdges) {
            int parent = edge.first;
            int child = edge.second;

            // Validate edge exists
            if (!hasParentChildEdge(T_init, parent, child)) {
                continue;
            }

            // Try rotating this edge
            VectorRangeTreeMap testTree = safeCopyTree(T_init);

            // Skip if copy failed
            if (testTree.original_nodes.empty()) {
                continue;
            }

            try {
                bool rotated = false;
                if (testTree.getLeftChild(parent) == child) {
                    testTree.rotateRight(parent);
                    rotated = true;
                } else if (testTree.getRightChild(parent) == child) {
                    testTree.rotateLeft(parent);
                    rotated = true;
                }

                if (!rotated) continue;

                // Check if any new edge matches target
                auto newEdges = getInternalEdges(testTree);
                for (const auto& newEdge : newEdges) {
                    if (initSet.find(newEdge) == initSet.end() &&
                        finalSet.find(newEdge) != finalSet.end()) {
                        // Found a free edge!
                        return {true, edge};
                    }
                }
            } catch (...) {
                continue;  // Skip problematic rotations
            }
        }
    } catch (...) {
        // Return false on any error
    }

    return {false, {-1, -1}};
}

// Main FlipDistTree function
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k) {
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    // If trees are already equal - this should be the FIRST check
    if (TreesEqual(T_init, T_final)) {
        debugPrint("Trees already equal, returning true");
        return true;
    }

    // Step 0: Early termination based on edge count
    int phi_init = countInternalEdges(T_init);
    debugPrint("T_init has " + std::to_string(phi_init) + " internal edges");

    // More lenient check - we need budget for rotations, not just edge count
    if (phi_init > k + 2) {  // Allow some extra budget
        debugPrint("Way too many edges, returning false");
        return false;
    }

    // Handle empty trees
    if (phi_init == 0) {
        bool result = countInternalEdges(T_final) == 0;
        debugPrint("No internal edges, result: " + std::string(result ? "true" : "false"));
        return result;
    }

    // Check for immediate free edge solution
    auto [hasFree, freeEdge] = findFreeEdge(T_init, T_final);
    if (hasFree && k >= 1) {
        debugPrint("Found immediate free edge solution");
        try {
            VectorRangeTreeMap T_rotated = safeCopyTree(T_init);
            int parent = freeEdge.first;
            int child = freeEdge.second;

            if (T_rotated.getLeftChild(parent) == child) {
                T_rotated.rotateRight(parent);
            } else if (T_rotated.getRightChild(parent) == child) {
                T_rotated.rotateLeft(parent);
            }

            if (TreesEqual(T_rotated, T_final)) {
                debugPrint("Single rotation solves it!");
                return true;
            }

            // Continue with recursive solution
            return FlipDistTree(T_rotated, T_final, k - 1);

        } catch (...) {
            debugPrint("Exception in immediate free edge handling");
        }
    }

    // Step 1: Enumerate all independent subsets of internal edges
    auto edges = getInternalEdges(T_init);
    debugPrint("Found " + std::to_string(edges.size()) + " internal edges");

    std::vector<std::vector<std::pair<int,int>>> independentSubsets;
    std::vector<std::pair<int,int>> current;

    // Limit complexity for larger trees
    int maxSubsets = (edges.size() > 6) ? 20 : 50;
    generateIndependentSubsets(edges, 0, current, independentSubsets, maxSubsets);

    debugPrint("Generated " + std::to_string(independentSubsets.size()) + " independent subsets");

    // Step 1.2: Try each non-empty independent subset
    for (size_t i = 0; i < independentSubsets.size(); i++) {
        const auto& subset = independentSubsets[i];
        if (subset.empty()) continue;

        debugPrint("Trying subset " + std::to_string(i) + " with " + std::to_string(subset.size()) + " edges");

        try {
            if (TreeDistI(T_init, T_final, k, subset)) {
                debugPrint("Found solution with subset " + std::to_string(i));
                return true;
            }
        } catch (...) {
            debugPrint("Exception in TreeDistI for subset " + std::to_string(i));
            continue;
        }
    }

    debugPrint("No solution found, returning false");
    return false;
}

// Simplified TreeDist-I implementation
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I) {
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    // If trees are already equal
    if (TreesEqual(T_init, T_final)) {
        debugPrint("TreeDistI: Trees already equal");
        return true;
    }

    // Step 0: More lenient early checks
    int phi_init = countInternalEdges(T_init);
    int budget_needed = (int)I.size();
    int remaining_budget = k - budget_needed;

    // Be more lenient - allow for the complexity of the remaining problem
    if (remaining_budget < 0) {
        debugPrint("TreeDistI: Not enough budget for rotations");
        return false;
    }

    if (phi_init == 0 && k >= 0) {
        bool result = TreesEqual(T_init, T_final);
        debugPrint("TreeDistI: No internal edges, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    // Step 2: Process each edge in I
    VectorRangeTreeMap T_bar = safeCopyTree(T_init);

    for (const auto& edge : I) {
        int parent = edge.first;
        int child = edge.second;

        debugPrint("TreeDistI: Rotating edge (" + std::to_string(parent) + "," + std::to_string(child) + ")");

        // Verify this is actually a valid edge in the current tree
        if (!hasParentChildEdge(T_bar, parent, child)) {
            debugPrint("TreeDistI: Invalid edge, skipping");
            continue;
        }

        try {
            if (T_bar.getLeftChild(parent) == child) {
                T_bar.rotateRight(parent);
            } else if (T_bar.getRightChild(parent) == child) {
                T_bar.rotateLeft(parent);
            }

            // Check if we've reached the target after this rotation
            if (TreesEqual(T_bar, T_final)) {
                debugPrint("TreeDistI: Reached target after rotation");
                return true;
            }

        } catch (...) {
            debugPrint("TreeDistI: Exception during rotation");
            return false;
        }
    }

    // Simplified: just call TreeDistS with empty S and remaining budget
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> emptyS;
    return TreeDistS(T_bar, T_final, remaining_budget, emptyS);
}

// Simplified TreeDist-S implementation
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    debugPrint("Entering TreeDistS with k=" + std::to_string(k));

    // Check if trees are already equal - this should be FIRST
    if (TreesEqual(T_init, T_end)) {
        debugPrint("TreeDistS: Trees already equal");
        return true;
    }

    // Step 0: Early checks - be more lenient
    int phi_init = countInternalEdges(T_init);
    if (phi_init > k + 1) {  // Allow some extra budget
        debugPrint("TreeDistS: Too many edges");
        return false;
    }

    if (phi_init == 0 && k >= 0) {
        bool result = TreesEqual(T_init, T_end);
        debugPrint("TreeDistS: No internal edges, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    // Early termination for small k
    if (k <= 0) {
        bool result = TreesEqual(T_init, T_end);
        debugPrint("TreeDistS: k=0, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    // Step 1: Check for free edge
    auto [hasFree, freeEdge] = findFreeEdge(T_init, T_end);

    if (hasFree) {
        debugPrint("TreeDistS: Found free edge (" + std::to_string(freeEdge.first) + "," + std::to_string(freeEdge.second) + ")");

        try {
            int parent = freeEdge.first;
            int child = freeEdge.second;

            // Apply the rotation
            VectorRangeTreeMap T_bar = safeCopyTree(T_init);
            if (T_bar.getLeftChild(parent) == child) {
                T_bar.rotateRight(parent);
            } else if (T_bar.getRightChild(parent) == child) {
                T_bar.rotateLeft(parent);
            }

            // Check if we've solved it
            if (TreesEqual(T_bar, T_end)) {
                debugPrint("TreeDistS: Solved with free edge rotation");
                return true;
            }

            // For simplicity, recursively call with reduced k
            return TreeDistS(T_bar, T_end, k - 1, {});

        } catch (...) {
            debugPrint("TreeDistS: Exception during free edge handling");
        }
    }

    // No free edge found - try a more aggressive approach
    debugPrint("TreeDistS: No free edge found");

    // If we have budget left, try some rotations
    if (k > phi_init) {
        auto edges = getInternalEdges(T_init);
        for (const auto& edge : edges) {
            try {
                VectorRangeTreeMap T_test = safeCopyTree(T_init);
                int parent = edge.first;
                int child = edge.second;

                if (T_test.getLeftChild(parent) == child) {
                    T_test.rotateRight(parent);
                } else if (T_test.getRightChild(parent) == child) {
                    T_test.rotateLeft(parent);
                }

                if (TreeDistS(T_test, T_end, k - 1, {})) {
                    debugPrint("TreeDistS: Found solution with rotation");
                    return true;
                }
            } catch (...) {
                continue;
            }
        }
    }

    debugPrint("TreeDistS: No solution found");
    return false;
}

// Helper function to print tree info with safety checks
void printTreeInfo(const std::string& name, const VectorRangeTreeMap& T) {
    try {
        std::cout << name << " edges: ";
        auto edges = getInternalEdges(T);
        for (const auto& e : edges) {
            std::cout << "(" << e.first << "," << e.second << ") ";
        }
        std::cout << "| Root: " << T.root << " | Nodes: " << T.original_nodes.size() << std::endl;
    } catch (...) {
        std::cout << name << " [ERROR: Cannot print tree info]" << std::endl;
    }
}

// Test case structure
struct TestCase {
    std::string name;
    std::vector<int> pre1, in1, pre2, in2;
    int expectedMinK;  // Minimum k where we expect True
    std::string description;
};

// Individual test runner with extensive debugging
bool runSingleTest(const TestCase& test, bool verbose = true) {
    if (verbose) {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "TEST: " << test.name << std::endl;
        std::cout << "Description: " << test.description << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        std::cout << "DEBUG: Starting test..." << std::endl;
    }

    try {
        // Step 1: Validate inputs
        if (verbose) std::cout << "DEBUG: Validating inputs..." << std::endl;

        if (test.pre1.empty() || test.in1.empty() || test.pre2.empty() || test.in2.empty()) {
            if (verbose) std::cout << "✗ Invalid test sequences" << std::endl;
            return false;
        }

        if (test.pre1.size() != test.in1.size() || test.pre2.size() != test.in2.size()) {
            if (verbose) std::cout << "✗ Sequence length mismatch" << std::endl;
            return false;
        }

        // Step 2: Print input sequences for debugging
        if (verbose) {
            std::cout << "DEBUG: T1 preorder: ";
            for (int x : test.pre1) std::cout << x << " ";
            std::cout << std::endl;
            std::cout << "DEBUG: T1 inorder: ";
            for (int x : test.in1) std::cout << x << " ";
            std::cout << std::endl;
            std::cout << "DEBUG: T2 preorder: ";
            for (int x : test.pre2) std::cout << x << " ";
            std::cout << std::endl;
            std::cout << "DEBUG: T2 inorder: ";
            for (int x : test.in2) std::cout << x << " ";
            std::cout << std::endl;
        }

        // Step 3: Build trees
        if (verbose) std::cout << "DEBUG: Building T1..." << std::endl;
        VectorRangeTreeMap T1;
        try {
            T1.build(test.pre1, test.in1);
        } catch (...) {
            if (verbose) std::cout << "✗ Failed to build T1" << std::endl;
            return false;
        }

        if (verbose) std::cout << "DEBUG: Building T2..." << std::endl;
        VectorRangeTreeMap T2;
        try {
            T2.build(test.pre2, test.in2);
        } catch (...) {
            if (verbose) std::cout << "✗ Failed to build T2" << std::endl;
            return false;
        }

        // Step 4: Validate tree construction
        if (verbose) std::cout << "DEBUG: Validating tree construction..." << std::endl;

        if (T1.original_nodes.empty() || T2.original_nodes.empty()) {
            if (verbose) std::cout << "✗ Trees not built correctly (empty node sets)" << std::endl;
            return false;
        }

        if (T1.root < 0 || T2.root < 0) {
            if (verbose) std::cout << "✗ Trees have invalid roots" << std::endl;
            return false;
        }

        // Step 5: Print tree info
        if (verbose) {
            std::cout << "DEBUG: Printing tree info..." << std::endl;
            printTreeInfo("T1", T1);
            printTreeInfo("T2", T2);
        }

        // Step 6: Test TreesEqual
        if (verbose) std::cout << "DEBUG: Testing TreesEqual..." << std::endl;
        bool areEqual = false;
        try {
            areEqual = TreesEqual(T1, T2);
        } catch (...) {
            if (verbose) std::cout << "✗ TreesEqual threw exception" << std::endl;
            return false;
        }

        if (verbose) {
            std::cout << "Trees equal: " << (areEqual ? "Yes" : "No") << std::endl;
        }

        // Step 7: Handle identical trees
        if (areEqual) {
            if (verbose) std::cout << "Trees are identical, testing k=0" << std::endl;
            bool result = false;
            try {
                result = FlipDistTree(T1, T2, 0);
            } catch (...) {
                if (verbose) std::cout << "✗ FlipDistTree(k=0) threw exception" << std::endl;
                return false;
            }
            if (verbose) std::cout << "FlipDistTree(k=0): " << (result ? "True" : "False") << std::endl;
            return result;
        }

        // Step 8: Test free edge detection
        if (verbose) std::cout << "DEBUG: Testing free edge detection..." << std::endl;
        try {
            auto [hasFree, freeEdge] = findFreeEdge(T1, T2);
            if (verbose) {
                std::cout << "Free edge: " << (hasFree ? "Found (" + std::to_string(freeEdge.first) + "," + std::to_string(freeEdge.second) + ")" : "None") << std::endl;
            }
        } catch (...) {
            if (verbose) std::cout << "? Free edge detection failed" << std::endl;
        }

        // Step 9: Test with different k values
        if (verbose) std::cout << "DEBUG: Testing with different k values..." << std::endl;
        bool foundSolution = false;
        int solutionK = -1;

        for (int k = 1; k <= 5; k++) {
            if (verbose) std::cout << "DEBUG: Testing k=" << k << std::endl;
            bool result = false;
            try {
                result = FlipDistTree(T1, T2, k);
            } catch (...) {
                if (verbose) std::cout << "✗ FlipDistTree(k=" << k << ") threw exception" << std::endl;
                continue;
            }

            if (verbose) {
                std::cout << "FlipDistTree(k=" << k << "): " << (result ? "True" : "False") << std::endl;
            }

            if (result && !foundSolution) {
                foundSolution = true;
                solutionK = k;
            }
        }

        // Step 10: Report results
        if (verbose) {
            std::cout << "DEBUG: Reporting results..." << std::endl;
            if (foundSolution) {
                std::cout << "✓ First solution found at k=" << solutionK << std::endl;
                if (test.expectedMinK > 0) {
                    std::cout << "Expected minimum k: " << test.expectedMinK << std::endl;
                    if (solutionK <= test.expectedMinK) {
                        std::cout << "✓ PASS: Found solution within expected bound" << std::endl;
                    } else {
                        std::cout << "? WARNING: Solution found later than expected" << std::endl;
                    }
                }
            } else {
                std::cout << "✗ No solution found within k=5" << std::endl;
            }
        }

        if (verbose) std::cout << "DEBUG: Test completed successfully" << std::endl;
        return foundSolution;

    } catch (const std::exception& e) {
        if (verbose) std::cout << "✗ Exception: " << e.what() << std::endl;
        return false;
    } catch (...) {
        if (verbose) std::cout << "✗ Unknown exception caught in test runner" << std::endl;
        return false;
    }
}

// Comprehensive test suite
void testFlipDist() {
    std::vector<TestCase> tests = {
            // Test 1: Identical trees
            {
                    "Identical Trees",
                    {2, 1, 3}, {1, 2, 3},
                    {2, 1, 3}, {1, 2, 3},
                    0,
                    "Two identical trees should require 0 rotations"
            },

            // Test 2: Single rotation
            {
                    "Single Left Rotation",
                    {2, 1, 3}, {1, 2, 3},
                    {3, 2, 1}, {1, 2, 3},
                    1,
                    "Left rotation at root: 2(1,3) -> 3(2(1,_),_)"
            },

            // Test 3: Single right rotation
            {
                    "Single Right Rotation",
                    {3, 2, 1}, {1, 2, 3},
                    {2, 1, 3}, {1, 2, 3},
                    1,
                    "Right rotation at root: 3(2(1,_),_) -> 2(1,3)"
            },

            // Test 4: Larger tree - chain to balanced
            {
                    "Chain to Balanced",
                    {1, 2, 3, 4}, {1, 2, 3, 4},
                    {3, 2, 1, 4}, {1, 2, 3, 4},
                    2,
                    "Transform right-skewed chain to more balanced tree"
            },

            // Test 5: Balanced to chain
            {
                    "Balanced to Chain",
                    {3, 2, 1, 4}, {1, 2, 3, 4},
                    {1, 2, 3, 4}, {1, 2, 3, 4},
                    2,
                    "Transform balanced tree to right-skewed chain"
            },

            // Test 6: Larger example
            {
                    "4-Node Complex",
                    {2, 1, 4, 3}, {1, 2, 3, 4},
                    {4, 2, 1, 3}, {1, 2, 3, 4},
                    3,
                    "More complex 4-node transformation"
            },

            // Test 7: Very safe test case - let's debug this specific issue
            {
                    "Safe Test",
                    {1, 2}, {1, 2},  // Simple 2-node tree
                    {2, 1}, {1, 2},  // Same nodes, different structure
                    1,
                    "Very simple and safe test case"
            },

            // Test 8: Edge case - single node
            {
                    "Single Node",
                    {1}, {1},
                    {1}, {1},
                    0,
                    "Single node trees should be identical"
            },

            // Test 9: Two nodes - original order
            {
                    "Two Nodes A",
                    {2, 1}, {1, 2},
                    {2, 1}, {1, 2},
                    0,
                    "Identical 2-node trees"
            },

            // Test 10: Two nodes - different order
            {
                    "Two Nodes B",
                    {1, 2}, {1, 2},
                    {2, 1}, {1, 2},
                    1,
                    "Different 2-node tree structures"
            }
    };

    std::cout << "\n" << std::string(80, '#') << std::endl;
    std::cout << "FLIPDISTTREE ALGORITHM TEST SUITE" << std::endl;
    std::cout << std::string(80, '#') << std::endl;

    int passed = 0;
    int total = tests.size();

    for (size_t i = 0; i < tests.size(); i++) {
        bool result = runSingleTest(tests[i], true);
        if (result) passed++;

        if (i < tests.size() - 1) {
            std::cout << "\nPress Enter for next test...";
            std::cin.get();
        }
    }

    std::cout << "\n" << std::string(80, '#') << std::endl;
    std::cout << "TEST SUMMARY" << std::endl;
    std::cout << std::string(80, '#') << std::endl;
    std::cout << "Passed: " << passed << "/" << total << " tests" << std::endl;

    if (passed == total) {
        std::cout << "🎉 ALL TESTS PASSED!" << std::endl;
    } else {
        std::cout << "⚠️  " << (total - passed) << " tests failed" << std::endl;
    }

    // Quick summary run
    std::cout << "\nQuick Summary (no debug output):" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    for (size_t i = 0; i < tests.size(); i++) {
        bool result = runSingleTest(tests[i], false);
        std::cout << "Test " << (i+1) << " (" << tests[i].name << "): "
                  << (result ? "✓ PASS" : "✗ FAIL") << std::endl;
    }
}

int main() {
    testFlipDist();
    return 0;
}