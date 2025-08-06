#include "A_tree.h"
#include <vector>
#include <iostream>
#include <utility>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <functional>

struct Edge {
    int u, v;
    bool operator==(const Edge &o) const { return u == o.u && v == o.v; }
};

namespace std {
    template<>
    struct hash<Edge> {
        size_t operator()(const Edge &e) const {
            return hash<long long>()(((long long)e.u << 32) ^ (unsigned int)e.v);
        }
    };
}

// Forward declarations
bool TreeDistI(VectorRangeTreeMap T1, const VectorRangeTreeMap& T2, int k, const std::vector<Edge>& I);
bool TreeDistS(VectorRangeTreeMap T1, const VectorRangeTreeMap& T2, int k, std::vector<std::pair<Edge, Edge>> S);

// Count internal nodes (non-leaf nodes) - this is the correct φ(T)
int countInternalNodes(const VectorRangeTreeMap &T) {
    int count = 0;
    for (int node : T.original_nodes) {
        bool isLeaf = (T.getLeftChild(node) == VectorRangeTreeMap::NO_CHILD &&
                       T.getRightChild(node) == VectorRangeTreeMap::NO_CHILD);
        if (!isLeaf) {
            count++;
        }
    }
    return count;
}

// Get all edges where parent is internal (non-leaf)
std::vector<Edge> getInternalEdges(const VectorRangeTreeMap &T) {
    std::unordered_set<std::pair<int, int>, PairHash, PairEq> all;
    T.collectEdges(T.root, all);

    std::vector<Edge> result;
    for (auto &[u, v] : all) {
        // Include edge if parent (u) is internal (has at least one child)
        bool parent_is_internal = (T.getLeftChild(u) != VectorRangeTreeMap::NO_CHILD ||
                                   T.getRightChild(u) != VectorRangeTreeMap::NO_CHILD);
        if (parent_is_internal) {
            result.push_back({u, v});
        }
    }
    return result;
}

// Get rotatable edges (edges where we can actually perform rotations)
std::vector<Edge> getRotatableEdges(const VectorRangeTreeMap &T) {
    return getInternalEdges(T); // Same as internal edges for now
}

// Check if a set of edges is independent (no shared endpoints)
bool isIndependent(const std::vector<Edge> &edges) {
    std::set<int> used;
    for (auto &e : edges) {
        if (used.count(e.u) || used.count(e.v)) return false;
        used.insert(e.u);
        used.insert(e.v);
    }
    return true;
}

// Generate all independent subsets of edges
std::vector<std::vector<Edge>> enumerateIndependentSubsets(const std::vector<Edge> &edges) {
    std::vector<std::vector<Edge>> result;
    int n = edges.size();

    // Include empty set
    result.push_back({});

    // Include non-empty subsets
    for (int mask = 1; mask < (1 << n); ++mask) {
        std::vector<Edge> subset;
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) subset.push_back(edges[i]);
        }
        if (isIndependent(subset)) result.push_back(subset);
    }
    return result;
}

// Debug function
void debugTree(const VectorRangeTreeMap& T, const std::string& name) {
    std::cout << "\n=== DEBUG: " << name << " ===" << std::endl;
    std::cout << "Root: " << T.root << std::endl;

    // Print original nodes
    std::vector<int> nodes(T.original_nodes.begin(), T.original_nodes.end());
    std::sort(nodes.begin(), nodes.end());

    std::cout << "Original nodes: ";
    for (int n : nodes) std::cout << n << " ";
    std::cout << std::endl;

    // Print structure
    for (int n : nodes) {
        auto range = T.getRange(n);
        int left = T.getLeftChild(n);
        int right = T.getRightChild(n);
        int parent = T.getParent(n);

        std::cout << "Node " << n << ": range(" << range.first << "," << range.second
                  << ") left=" << left << " right=" << right << " parent=" << parent;

        // Check if leaf
        bool isLeaf = (left == VectorRangeTreeMap::NO_CHILD && right == VectorRangeTreeMap::NO_CHILD);
        if (isLeaf) std::cout << " [LEAF]";
        std::cout << std::endl;
    }

    int internal_count = countInternalNodes(T);
    std::cout << "φ(T) = " << internal_count << " internal nodes" << std::endl;

    auto rotatable_edges = getRotatableEdges(T);
    std::cout << "Rotatable edges (" << rotatable_edges.size() << "): ";
    for (auto& e : rotatable_edges) {
        std::cout << "(" << e.u << "," << e.v << ") ";
    }
    std::cout << std::endl;
}

bool FlipDistTree(const VectorRangeTreeMap& T1, const VectorRangeTreeMap& T2, int k) {
    std::cout << "\n[FlipDistTree] Starting with k=" << k << std::endl;

    // Step 0: If trees are equal, return true
    if (TreesEqual(T1, T2)) {
        std::cout << "[FlipDistTree] ✓ Trees are equal, returning TRUE" << std::endl;
        return true;
    }

    // φ(T1) = number of internal nodes
    int phi_T1 = countInternalNodes(T1);
    std::cout << "[FlipDistTree] φ(T1) = " << phi_T1 << " internal nodes" << std::endl;

    if (phi_T1 > k) {
        std::cout << "[FlipDistTree] φ(T1) > k, returning FALSE" << std::endl;
        return false;
    }

    // Step 1: Enumerate all independent subsets of rotatable edges
    auto rotatableEdges = getRotatableEdges(T1);
    std::cout << "[FlipDistTree] Rotatable edges: " << rotatableEdges.size() << std::endl;

    auto subsets = enumerateIndependentSubsets(rotatableEdges);
    std::cout << "[FlipDistTree] Trying " << subsets.size() << " independent subsets" << std::endl;

    // Step 1.2: For each non-empty subset I
    for (auto& I : subsets) {
        if (I.empty()) continue;

        std::cout << "\n[FlipDistTree] → Trying subset of size " << I.size() << ": ";
        for (auto& e : I) std::cout << "(" << e.u << "," << e.v << ") ";
        std::cout << std::endl;

        // Make a copy since TreeDistI modifies the tree
        VectorRangeTreeMap T1_copy = T1;

        if (TreeDistI(T1_copy, T2, k, I)) {
            std::cout << "[FlipDistTree] ✓✓✓ FOUND SOLUTION!" << std::endl;
            return true;
        }
    }

    std::cout << "[FlipDistTree] No solution found" << std::endl;
    return false;
}

bool TreeDistI(VectorRangeTreeMap T1, const VectorRangeTreeMap& T2, int k, const std::vector<Edge>& I) {
    std::cout << "\n[TreeDistI] called with k=" << k << ", |I|=" << I.size() << std::endl;

    // CRITICAL: Check if we're already equal before any processing
    if (TreesEqual(T1, T2)) {
        std::cout << "[TreeDistI] ✓✓✓ Trees already equal! SUCCESS!" << std::endl;
        return true;
    }

    // Step 0: Check preconditions
    int phi_T1 = countInternalNodes(T1);
    int remaining_k = k - (int)I.size();

    std::cout << "[TreeDistI] φ(T1) = " << phi_T1 << " internal nodes, k - |I| = " << remaining_k << std::endl;

    // Allow the case where phi_T1 == remaining_k, since we can still perform exactly that many rotations
    if (phi_T1 > remaining_k) {
        std::cout << "[TreeDistI] φ(T1) > k - |I|, returning FALSE" << std::endl;
        return false;
    }

    if (phi_T1 == 0 && remaining_k >= 0) {
        std::cout << "[TreeDistI] φ(T1) = 0, checking equality" << std::endl;
        bool equal = TreesEqual(T1, T2);
        std::cout << "[TreeDistI] Trees equal? " << (equal ? "TRUE" : "FALSE") << std::endl;
        return equal;
    }

    // Step 2: For each edge e ∈ I, rotate e
    for (const auto& e : I) {
        int u = e.u, v = e.v; // parent, child
        std::cout << "\n[TreeDistI] → Processing edge (" << u << "," << v << ")" << std::endl;

        // Step 2.1: Perform rotation based on which child v is
        bool rotated = false;
        if (T1.getLeftChild(u) == v) {
            std::cout << "[TreeDistI]   Rotating RIGHT at " << u << " (pulling up left child " << v << ")" << std::endl;
            T1.rotateRight(u);
            rotated = true;
        } else if (T1.getRightChild(u) == v) {
            std::cout << "[TreeDistI]   Rotating LEFT at " << u << " (pulling up right child " << v << ")" << std::endl;
            T1.rotateLeft(u);
            rotated = true;
        }

        if (!rotated) {
            std::cout << "[TreeDistI] ERROR: Could not rotate edge (" << u << "," << v << ")" << std::endl;
            return false;
        }

        // CRITICAL: Check if we're done immediately after each rotation
        if (TreesEqual(T1, T2)) {
            std::cout << "[TreeDistI] ✓✓✓ TREES EQUAL after rotation! SUCCESS!" << std::endl;
            return true;
        }
    }

    //just check if we're equal after all rotations
    if (TreesEqual(T1, T2)) {
        std::cout << "[TreeDistI] ✓✓✓ TREES EQUAL after all rotations! SUCCESS!" << std::endl;
        return true;
    }

    std::cout << "[TreeDistI] Trees not equal after rotations, would call TreeDistS" << std::endl;
    return false;
}

bool TreeDistS(VectorRangeTreeMap T1, const VectorRangeTreeMap& T2, int k, std::vector<std::pair<Edge, Edge>> S) {
    // Simplified implementation for now
    return TreesEqual(T1, T2);
}

// Test the correct single rotation
void testCorrectRotation() {
    std::cout << "\n=== TESTING CORRECT ROTATION ===" << std::endl;

    VectorRangeTreeMap T1, T2;

    // T1: Tree with structure 3(1,2) - preorder {3,1,2}, inorder {1,3,2}
    T1.build({3, 1, 2}, {1, 3, 2});

    // T2: Tree with structure 1(NULL,3(2,NULL)) - preorder {1,3,2}, inorder {1,2,3}
    T2.build({1, 3, 2}, {1, 2, 3});

    debugTree(T1, "T1 (source)");
    debugTree(T2, "T2 (target)");

    // The correct rotation should be: rotate RIGHT at node 3
    // This pulls up node 1, making it the new root
    std::cout << "\n--- Manual rotation test ---" << std::endl;
    VectorRangeTreeMap T1_manual = T1;
    T1_manual.rotateRight(3);

    debugTree(T1_manual, "T1 after manual RIGHT rotation at 3");

    std::cout << "Trees equal after manual rotation? " << (TreesEqual(T1_manual, T2) ? "YES" : "NO") << std::endl;

    // Compare serializations
    std::cout << "\nSerialization comparison:" << std::endl;
    std::cout << "T1 manual: " << treeToString(T1_manual) << std::endl;
    std::cout << "T2 target: " << treeToString(T2) << std::endl;

    // Now test the algorithm
    std::cout << "\n--- Algorithm test ---" << std::endl;
    bool result = FlipDistTree(T1, T2, 1);
    std::cout << "Algorithm result (k=1): " << (result ? "SUCCESS" : "FAIL") << std::endl;
}

int main() {
    std::cout << "=== Running Fixed FlipDistTree Tests ===" << std::endl;

    // Test the correct rotation first
    testCorrectRotation();

    std::cout << "\n" << std::string(60, '=') << std::endl;

    VectorRangeTreeMap T1, T2;

    // Test 1: Identical trees
    std::cout << "\n=== Test 1: Identical Trees ===" << std::endl;
    T1.build({3, 1, 2}, {1, 3, 2});
    T2.build({3, 1, 2}, {1, 3, 2});
    std::cout << "Test 1 (identical): " << (FlipDistTree(T1, T2, 0) ? "✓ PASS" : "✗ FAIL") << std::endl;

    // Test 2: One rotation needed
    std::cout << "\n=== Test 2: One Rotation ===" << std::endl;
    T1.build({3, 1, 2}, {1, 3, 2});
    T2.build({1, 3, 2}, {1, 2, 3});

    debugTree(T1, "T1 (source)");
    debugTree(T2, "T2 (target)");

    // Test with k=1 (should work now)
    std::cout << "\n--- Testing with k=1 ---" << std::endl;
    bool result1 = FlipDistTree(T1, T2, 1);
    std::cout << "Test 2 (k=1): " << (result1 ? "✓ PASS" : "✗ FAIL") << std::endl;

    return 0;
}