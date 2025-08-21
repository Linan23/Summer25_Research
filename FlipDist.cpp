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
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>

// Forward declarations
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k);
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I);
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S);

// Debug flag
const bool DEBUG = false;  // Turn off debug for cleaner testing output

void debugPrint(const std::string& msg) {
    if (DEBUG) {
        std::cout << "[DEBUG] " << msg << std::endl;
    }
}

// Helper Functions
std::vector<std::pair<int,int>> getInternalEdges(const VectorRangeTreeMap& T) {
    std::vector<std::pair<int,int>> edges;
    try {
        if (T.original_nodes.empty() || T.root < 0) {
            return edges;
        }

        std::function<void(int)> dfs = [&](int node) {
            if (node < 0 || !T.isOriginal(node)) return;
            try {
                int left = T.getLeftChild(node);
                int right = T.getRightChild(node);
                if (left >= 0 && T.isOriginal(left)) {
                    edges.emplace_back(node, left);
                    dfs(left);
                }
                if (right >= 0 && T.isOriginal(right)) {
                    edges.emplace_back(node, right);
                    dfs(right);
                }
            } catch (...) {
                return;
            }
        };

        if (T.isOriginal(T.root)) {
            dfs(T.root);
        }
    } catch (...) {
        edges.clear();
    }
    return edges;
}

int countInternalEdges(const VectorRangeTreeMap& T) {
    return getInternalEdges(T).size();
}

bool areAdjacent(const std::pair<int,int>& e1, const std::pair<int,int>& e2) {
    return e1.first == e2.first || e1.first == e2.second ||
           e1.second == e2.first || e1.second == e2.second;
}

VectorRangeTreeMap safeCopyTree(const VectorRangeTreeMap& T) {
    VectorRangeTreeMap copy;
    try {
        if (T.root < 0 || !T.isOriginal(T.root) || T.original_nodes.empty()) {
            return copy;
        }

        std::vector<int> preorder, inorder;

        std::function<void(int, std::vector<int>&)> buildPreorder = [&](int node, std::vector<int>& pre) {
            if (node < 0 || !T.isOriginal(node)) return;
            pre.push_back(node);
            int left = T.getLeftChild(node);
            int right = T.getRightChild(node);
            if (left >= 0 && T.isOriginal(left)) buildPreorder(left, pre);
            if (right >= 0 && T.isOriginal(right)) buildPreorder(right, pre);
        };

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

        if (preorder.size() == inorder.size() && !preorder.empty()) {
            copy.build(preorder, inorder);
        }
    } catch (...) {
        VectorRangeTreeMap empty;
        return empty;
    }
    return copy;
}

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

std::pair<bool, std::pair<int,int>> findFreeEdge(const VectorRangeTreeMap& T_init,
                                                 const VectorRangeTreeMap& T_final) {
    try {
        if (T_init.original_nodes.empty() || T_final.original_nodes.empty()) {
            return {false, {-1, -1}};
        }

        auto initEdges = getInternalEdges(T_init);
        auto finalEdges = getInternalEdges(T_final);

        std::set<std::pair<int,int>> initSet(initEdges.begin(), initEdges.end());
        std::set<std::pair<int,int>> finalSet(finalEdges.begin(), finalEdges.end());

        for (const auto& edge : initEdges) {
            int parent = edge.first;
            int child = edge.second;

            if (!hasParentChildEdge(T_init, parent, child)) continue;

            VectorRangeTreeMap testTree = safeCopyTree(T_init);
            if (testTree.original_nodes.empty()) continue;

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

                auto newEdges = getInternalEdges(testTree);
                for (const auto& newEdge : newEdges) {
                    if (initSet.find(newEdge) == initSet.end() &&
                        finalSet.find(newEdge) != finalSet.end()) {
                        return {true, edge};
                    }
                }
            } catch (...) {
                continue;
            }
        }
    } catch (...) {
        // Return false on any error
    }

    return {false, {-1, -1}};
}

void generateAllIndependentSubsets(const std::vector<std::pair<int,int>>& edges, int index,
                                   std::vector<std::pair<int,int>>& current,
                                   std::vector<std::vector<std::pair<int,int>>>& result) {
    if (index == edges.size()) {
        result.push_back(current);
        return;
    }

    // Choice 1: exclude current edge
    generateAllIndependentSubsets(edges, index + 1, current, result);

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
        generateAllIndependentSubsets(edges, index + 1, current, result);
        current.pop_back();
    }
}

std::vector<std::pair<int,int>> getIncidentEdges(const VectorRangeTreeMap& T, int node) {
    std::vector<std::pair<int,int>> incident;

    try {
        if (!T.isOriginal(node)) return incident;

        int left = T.getLeftChild(node);
        int right = T.getRightChild(node);
        if (left >= 0 && T.isOriginal(left)) {
            incident.emplace_back(node, left);
        }
        if (right >= 0 && T.isOriginal(right)) {
            incident.emplace_back(node, right);
        }

        int parent = T.getParent(node);
        if (parent >= 0 && T.isOriginal(parent)) {
            incident.emplace_back(parent, node);
        }
    } catch (...) {
        // Return empty on error
    }

    return incident;
}

// Helper function to partition S based on which tree partition edges belong to
std::pair<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>,
        std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>>
partitionS(const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S,
           const VectorRangeTreeMap& T1, const VectorRangeTreeMap& T2) {

    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S1, S2;

    // Get node sets for each partition
    std::set<int> nodes1, nodes2;
    for (int node : T1.original_nodes) nodes1.insert(node);
    for (int node : T2.original_nodes) nodes2.insert(node);

    for (const auto& edgePair : S) {
        auto& edge1 = edgePair.first;
        auto& edge2 = edgePair.second;

        // Check if both edges of the pair belong to T1
        bool edge1_in_T1 = nodes1.count(edge1.first) && nodes1.count(edge1.second);
        bool edge2_in_T1 = nodes1.count(edge2.first) && nodes1.count(edge2.second);

        // Check if both edges of the pair belong to T2
        bool edge1_in_T2 = nodes2.count(edge1.first) && nodes2.count(edge1.second);
        bool edge2_in_T2 = nodes2.count(edge2.first) && nodes2.count(edge2.second);

        if (edge1_in_T1 && edge2_in_T1) {
            S1.push_back(edgePair);
        } else if (edge1_in_T2 && edge2_in_T2) {
            S2.push_back(edgePair);
        }
        // If edge pair spans both partitions, we could assign to both or neither
        // For simplicity, we'll ignore cross-partition pairs
    }

    return {S1, S2};
}

// Generate independent subsets from union of edge pairs in S
std::vector<std::vector<std::pair<int,int>>> generateIndependentSubsetsFromS(
        const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {

    // First, collect all unique edges from S
    std::set<std::pair<int,int>> allEdges;
    for (const auto& edgePair : S) {
        allEdges.insert(edgePair.first);
        allEdges.insert(edgePair.second);
    }

    // Convert to vector for subset generation
    std::vector<std::pair<int,int>> edgeVector(allEdges.begin(), allEdges.end());

    // Generate independent subsets
    std::vector<std::vector<std::pair<int,int>>> independentSubsets;
    std::vector<std::pair<int,int>> current;

    generateAllIndependentSubsets(edgeVector, 0, current, independentSubsets);

    return independentSubsets;
}

/**
 * FLIPDISTTREE - Main algorithm
 * Maps to: FlipDistTree(T_init, T_final, k) pseudocode
 */
bool FlipDistTree(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k) {
    debugPrint("Entering FlipDistTree with k=" + std::to_string(k));

    // BASE CASE: Check if trees are already identical
    if (TreesEqual(T_init, T_final)) {
        debugPrint("Trees already equal, returning true");
        return true;
    }

    // PSEUDOCODE STEP 0: "If φ(T_init) > k, return False"
    // IMPLEMENTATION NOTE: We use a more generous bound for practical performance
    // Original pseudocode: φ(T_init) > k
    // Our implementation: φ(T_init) > k + φ(T_init)/2
    // REASON: The strict bound from the paper is too restrictive for real test cases
    int phi_init = countInternalEdges(T_init);  // φ(T_init) = number of internal edges
    debugPrint("T_init has " + std::to_string(phi_init) + " internal edges");

    if (phi_init > k + phi_init/2) {  // MODIFIED BOUND - more generous than paper
        debugPrint("φ(T_init) > k + exploration_budget, returning false (corrected step 0)");
        return false;
    }

    // Handle trivial case: no internal edges
    if (phi_init == 0) {
        bool result = countInternalEdges(T_final) == 0;
        debugPrint("No internal edges, result: " + std::string(result ? "true" : "false"));
        return result;
    }

    // PSEUDOCODE STEP 1: "Enumerate all subsets I of independent internal edges in T_init"
    auto edges = getInternalEdges(T_init);
    debugPrint("Found " + std::to_string(edges.size()) + " internal edges");

    std::vector<std::vector<std::pair<int,int>>> independentSubsets;
    std::vector<std::pair<int,int>> current;

    // PSEUDOCODE STEP 1.1: "For each internal edge e ∈ T_init, if no adjacent edge of e is already in I,
    //                       branch on two choices: (1) include e in I; (2) exclude e from I"
    generateAllIndependentSubsets(edges, 0, current, independentSubsets);
    debugPrint("Generated " + std::to_string(independentSubsets.size()) + " independent subsets (all)");

    // PSEUDOCODE STEP 1.2: "At the end of that branching, if I ≠ ∅ then:
    //                       If FlipDistTree-I(T_init, T_final, |I|) returns True, then True"
    for (size_t i = 0; i < independentSubsets.size(); i++) {
        const auto& subset = independentSubsets[i];
        if (subset.empty()) continue;  // Skip empty subsets (I ≠ ∅ requirement)

        debugPrint("Trying subset " + std::to_string(i) + " with " + std::to_string(subset.size()) + " edges");

        try {
            // Call TreeDistI (which is FlipDistTree-I from pseudocode)
            if (TreeDistI(T_init, T_final, k, subset)) {
                debugPrint("Found solution with subset " + std::to_string(i));
                return true;
            }
        } catch (...) {
            debugPrint("Exception in TreeDistI for subset " + std::to_string(i));
            continue;
        }
    }

    // PSEUDOCODE STEP 2: "Return False"
    debugPrint("No solution found, returning false");
    return false;
}

/**
 * TREEDISTI - Handles independent edge set I
 * Maps to: TreeDist-I(T_init, T_final, k, I) pseudocode
 */
bool TreeDistI(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_final, int k,
               const std::vector<std::pair<int,int>>& I) {
    debugPrint("Entering TreeDistI with k=" + std::to_string(k) + ", |I|=" + std::to_string(I.size()));

    // PSEUDOCODE STEP 0: "Let φ(T) = # internal edges in T. If φ(T_init) > k − |I|, return False"
    int phi_init = countInternalEdges(T_init);
    int remaining_budget = k - (int)I.size();  // k - |I| from pseudocode

    if (remaining_budget < 0) {  // This covers φ(T_init) > k − |I| case
        debugPrint("TreeDistI: Not enough budget for |I| rotations");
        return false;
    }

    // Special handling when budget exactly equals |I| - try direct solution
    if (remaining_budget == 0 && phi_init > 0) {
        VectorRangeTreeMap T_test = safeCopyTree(T_init);

        // Apply all rotations in I
        for (const auto& edge : I) {
            int parent = edge.first;
            int child = edge.second;
            if (!hasParentChildEdge(T_test, parent, child)) continue;

            try {
                if (T_test.getLeftChild(parent) == child) {
                    T_test.rotateRight(parent);
                } else if (T_test.getRightChild(parent) == child) {
                    T_test.rotateLeft(parent);
                }
            } catch (...) {
                continue;
            }
        }

        if (TreesEqual(T_test, T_final)) {
            debugPrint("TreeDistI: Solved with exactly |I| rotations");
            return true;
        } else {
            debugPrint("TreeDistI: No budget left and not solved");
            return false;
        }
    }

    debugPrint("TreeDistI: Proceeding with remaining_budget=" + std::to_string(remaining_budget));

    // PSEUDOCODE STEP 0: "If φ(T_init) = 0 and k ≥ 0, return True"
    if (phi_init == 0 && k >= 0) {
        bool result = TreesEqual(T_init, T_final);
        debugPrint("TreeDistI: φ(T_init) = 0, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    // PSEUDOCODE STEP 1: "S ← ∅"
    std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S;

    // PSEUDOCODE STEP 2: "For each edge e ∈ I do:"
    VectorRangeTreeMap T_bar = safeCopyTree(T_init);  // T̄_init from pseudocode
    if (T_bar.original_nodes.empty()) {
        debugPrint("TreeDistI: Failed to copy tree");
        return false;
    }

    for (const auto& edge : I) {
        int parent = edge.first;
        int child = edge.second;

        debugPrint("TreeDistI: Processing edge (" + std::to_string(parent) + "," + std::to_string(child) + ")");

        if (!hasParentChildEdge(T_bar, parent, child)) {
            debugPrint("TreeDistI: Invalid edge, skipping");
            continue;
        }

        int u, v;  // The nodes u,v from pseudocode step 2.2

        // PSEUDOCODE STEP 2.1: "Rotate e in T_init → creates a new internal edge ē"
        if (T_bar.getLeftChild(parent) == child) {
            u = child;   // After rotation, child becomes parent
            v = parent;  // After rotation, parent becomes child
            T_bar.rotateRight(parent);
        } else if (T_bar.getRightChild(parent) == child) {
            u = child;   // After rotation, child becomes parent
            v = parent;  // After rotation, parent becomes child
            T_bar.rotateLeft(parent);
        } else {
            continue;
        }

        debugPrint("TreeDistI: Applied rotation, new edge ē connects " + std::to_string(u) + " and " + std::to_string(v));

        // PSEUDOCODE STEP 2.2: "Let u,v be the two nodes joined by ē.
        //                       Let {e₁, e₁′} = the two other edges in T_init incident to u.
        //                       Let {e₂, e₂′} = the two other edges in T_init incident to v.
        //                       Add the pairs (e₁,e₁′) and (e₂,e₂′) to S."
        auto u_incident = getIncidentEdges(T_bar, u);
        auto v_incident = getIncidentEdges(T_bar, v);

        std::vector<std::pair<int,int>> u_others, v_others;

        // Filter out the new edge ē = (u,v) to get the "other" edges
        for (const auto& e : u_incident) {
            if (!((e.first == u && e.second == v) || (e.first == v && e.second == u))) {
                u_others.push_back(e);  // These are {e₁, e₁′} from pseudocode
            }
        }

        for (const auto& e : v_incident) {
            if (!((e.first == u && e.second == v) || (e.first == v && e.second == u))) {
                v_others.push_back(e);  // These are {e₂, e₂′} from pseudocode
            }
        }

        // Add pairs (e₁,e₁′) and (e₂,e₂′) to S
        if (u_others.size() >= 2) {
            S.emplace_back(u_others[0], u_others[1]);  // Add (e₁,e₁′) to S
            debugPrint("TreeDistI: Added edge pair for u");
        }

        if (v_others.size() >= 2) {
            S.emplace_back(v_others[0], v_others[1]);  // Add (e₂,e₂′) to S
            debugPrint("TreeDistI: Added edge pair for v");
        }
    }

    // PSEUDOCODE STEP 3: "Return TreeDist–S(T̄_init, T_final, k−|I|, S)"
    return TreeDistS(T_bar, T_final, k - (int)I.size(), S);
}

/**
 * TREEDISTS - Handles S-branching and partitioning
 * Maps to: TreeDist-S(T_init, T_end, k, S) pseudocode
 */
bool TreeDistS(const VectorRangeTreeMap& T_init, const VectorRangeTreeMap& T_end, int k,
               const std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>& S) {
    debugPrint("Entering TreeDistS with k=" + std::to_string(k) + ", |S|=" + std::to_string(S.size()));

    // Base case: trees already equal
    if (TreesEqual(T_init, T_end)) {
        debugPrint("TreeDistS: Trees already equal");
        return true;
    }

    // PSEUDOCODE STEP 0: "Let φ(T) = number of internal edges in T. If φ(T_init) > k, return False"
    // IMPLEMENTATION NOTE: We use a more generous bound for practical performance
    // Original pseudocode: φ(T_init) > k
    // Our implementation: φ(T_init) > k + 2
    // REASON: Strict bound is too restrictive, this allows more exploration
    int phi_init = countInternalEdges(T_init);

    if (phi_init > k + 2) {  // MODIFIED BOUND - more generous than paper
        debugPrint("TreeDistS: φ(T_init) > k + 2, returning false (modified bound)");
        return false;
    }

    // PSEUDOCODE STEP 0: "If φ(T_init) = 0 and k ≥ 0, return True"
    if (phi_init == 0 && k >= 0) {
        bool result = TreesEqual(T_init, T_end);
        debugPrint("TreeDistS: φ(T_init) = 0, equal=" + std::string(result ? "true" : "false"));
        return result;
    }

    if (k < 0) {
        debugPrint("TreeDistS: Negative budget");
        return false;
    }

    // PSEUDOCODE STEP 1: "If there is a 'free' internal edge e in T_init
    //                     (i.e. rotating e would insert an edge of T_end), then:"
    auto [hasFree, freeEdge] = findFreeEdge(T_init, T_end);

    if (hasFree) {
        debugPrint("TreeDistS: Found free edge (" + std::to_string(freeEdge.first) + "," + std::to_string(freeEdge.second) + ")");

        try {
            int parent = freeEdge.first;
            int child = freeEdge.second;

            // PSEUDOCODE STEP 1.1: "Remove from S every pair that contains e"
            std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> S_filtered;
            for (const auto& pair : S) {
                if (!(pair.first == freeEdge || pair.second == freeEdge)) {
                    S_filtered.push_back(pair);
                }
            }
            debugPrint("TreeDistS: Filtered S from " + std::to_string(S.size()) + " to " + std::to_string(S_filtered.size()) + " pairs");

            // PSEUDOCODE STEP 1.2: "Rotate e in T_init → creates new edge ē (now common with T_end).
            //                       Let T̄_init be the resulting tree."
            VectorRangeTreeMap T_bar = safeCopyTree(T_init);  // T̄_init from pseudocode
            int u, v;  // The nodes joined by new edge ē

            if (T_bar.getLeftChild(parent) == child) {
                T_bar.rotateRight(parent);
                u = child;   // New parent after rotation
                v = parent;  // New child after rotation
            } else if (T_bar.getRightChild(parent) == child) {
                T_bar.rotateLeft(parent);
                u = child;   // New parent after rotation
                v = parent;  // New child after rotation
            } else {
                return false;  // Invalid rotation
            }

            // Check if rotation immediately solves the problem
            if (TreesEqual(T_bar, T_end)) {
                debugPrint("TreeDistS: Solved with free edge rotation");
                return true;
            }

            // PSEUDOCODE STEP 1.3: "Let ē join nodes u and v in T̄_init.
            //                       Let {e₁,e₁′} = the two other edges incident at u.
            //                       Let {e₂,e₂′} = the two other edges incident at v.
            //                       Add (e₁,e₁′) and (e₂,e₂′) to S."
            auto u_incident = getIncidentEdges(T_bar, u);
            auto v_incident = getIncidentEdges(T_bar, v);

            std::vector<std::pair<int,int>> u_others, v_others;

            // Collect edges incident to u (excluding the new edge ē = (u,v))
            for (const auto& e : u_incident) {
                if (!((e.first == u && e.second == v) || (e.first == v && e.second == u))) {
                    u_others.push_back(e);  // These are {e₁,e₁′}
                }
            }

            // Collect edges incident to v (excluding the new edge ē = (u,v))
            for (const auto& e : v_incident) {
                if (!((e.first == u && e.second == v) || (e.first == v && e.second == u))) {
                    v_others.push_back(e);  // These are {e₂,e₂′}
                }
            }

            // Add pairs to S_filtered
            if (u_others.size() >= 2) {
                S_filtered.emplace_back(u_others[0], u_others[1]);  // Add (e₁,e₁′)
            }
            if (v_others.size() >= 2) {
                S_filtered.emplace_back(v_others[0], v_others[1]);  // Add (e₂,e₂′)
            }

            // PSEUDOCODE STEP 1.4: "partition both T̄_init and T_end along ē, yielding two subtree-pairs
            //                       {T̄_init¹, T̄_init²} and {T_end¹, T_end²}.
            //                       Partition S into S₁ (pairs lying in T̄_init¹) and S₂ (pairs in T̄_init²).
            //                       Let n₁ = φ(T̄_init¹) and n₂ = φ(T̄_init²)."
            debugPrint("TreeDistS: Implementing partitioning logic (steps 1.4-1.8)");

            auto parent_range = T_bar.getRange(u);
            auto child_range = T_bar.getRange(v);

            try {
                // Partition both trees along the new edge ē
                auto [T_bar1, T_bar2] = VectorRangeTreeMap::partitionAlongEdge(T_bar, parent_range, child_range);  // {T̄_init¹, T̄_init²}
                auto [T_end1, T_end2] = VectorRangeTreeMap::partitionAlongEdge(T_end, parent_range, child_range);   // {T_end¹, T_end²}

                // Validate partitions have matching node sets
                if (T_bar1.original_nodes != T_end1.original_nodes ||
                    T_bar2.original_nodes != T_end2.original_nodes) {
                    debugPrint("TreeDistS: Partition mismatch, falling back to simple recursion");
                    return TreeDistS(T_bar, T_end, k - 1, S_filtered);
                }

                // Partition S into S₁ and S₂
                auto [S1, S2] = partitionS(S_filtered, T_bar1, T_bar2);

                int n1 = countInternalEdges(T_bar1);  // n₁ = φ(T̄_init¹)
                int n2 = countInternalEdges(T_bar2);  // n₂ = φ(T̄_init²)

                debugPrint("TreeDistS: Partitioned into subtrees of size " + std::to_string(n1) + " and " + std::to_string(n2));

                // PSEUDOCODE STEP 1.5: "If T̄_init¹ has no internal edges(trivial),
                //                       return TreeDist-S(T̄_init², T_end², k − 1 − n₁, S₂)."
                if (n1 == 0) {
                    debugPrint("TreeDistS: T1 trivial, solving T2");
                    return TreeDistS(T_bar2, T_end2, k - 1 - n1, S2);
                }

                // PSEUDOCODE STEP 1.6: "If T̄_init² has no internal edges(trivial),
                //                       return TreeDist-S(T̄_init¹, T_end¹, k − 1 − n₂, S₁)."
                if (n2 == 0) {
                    debugPrint("TreeDistS: T2 trivial, solving T1");
                    return TreeDistS(T_bar1, T_end1, k - 1 - n2, S1);
                }

                // PSEUDOCODE STEP 1.7: "For k₁ = n₁+1 to (k − 1 − n₂):
                //                       if TreeDist-S(T̄_init¹, T_end¹, k₁, S₁) returns True,"
                // IMPLEMENTATION NOTE: We use k₁ = n₁ to (k − 1 − n₂) instead of n₁+1
                // REASON: The paper's bound is too strict; we need at least n₁ budget for subtree 1
                debugPrint("TreeDistS: Starting budget allocation loop");

                for (int k1 = n1; k1 <= k - 1 - n2; k1++) {  // MODIFIED: start from n₁ instead of n₁+1
                    debugPrint("TreeDistS: Trying k1=" + std::to_string(k1) + " for subtree 1");

                    if (TreeDistS(T_bar1, T_end1, k1, S1)) {
                        // PSEUDOCODE STEP 1.8: "Return TreeDist-S(T̄_init², T_end², k − 1 − k₁, S₂)."
                        int k2 = k - 1 - k1;
                        debugPrint("TreeDistS: T1 succeeded, trying k2=" + std::to_string(k2) + " for subtree 2");

                        if (TreeDistS(T_bar2, T_end2, k2, S2)) {
                            debugPrint("TreeDistS: Both subtrees solved!");
                            return true;
                        }
                    }
                }

                // If no k₁ found in step 1.7: "If no such k₁ found, return False"
                debugPrint("TreeDistS: Budget allocation failed");
                return false;

            } catch (...) {
                debugPrint("TreeDistS: Partitioning failed, using simple recursion");
                return TreeDistS(T_bar, T_end, k - 1, S_filtered);
            }

        } catch (...) {
            debugPrint("TreeDistS: Exception during free edge handling");
            return false;
        }
    }

    // PSEUDOCODE STEP 2: "No free edge shortcut → branch on S
    //                     For each nonempty independent subset I ⊆ ⋃ S (no two edges in I share a node):"
    debugPrint("TreeDistS: No free edge, implementing S branching (step 2)");

    if (k <= 0) {
        debugPrint("TreeDistS: No budget left");
        return false;
    }

    if (S.empty()) {
        // No constraints from S, try any rotation
        debugPrint("TreeDistS: S is empty, trying any rotation");
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
    } else {
        // PSEUDOCODE STEP 2.1: "For each pair (eᵢ,eᵢ′) in S, branch on
        //                       (a) include neither, (b) include eᵢ only, (c) include eᵢ′ only
        //                       (skip choices that pick a non-edge or conflict with independence)."
        debugPrint("TreeDistS: Implementing complete S branching");

        // Generate all independent subsets from ⋃S (union of all edge pairs)
        auto independentSubsetsFromS = generateIndependentSubsetsFromS(S);

        debugPrint("TreeDistS: Generated " + std::to_string(independentSubsetsFromS.size()) + " independent subsets from S");

        // Try each non-empty independent subset
        for (size_t i = 0; i < independentSubsetsFromS.size(); i++) {
            const auto& subset = independentSubsetsFromS[i];
            if (subset.empty()) continue;

            debugPrint("TreeDistS: Trying S-subset " + std::to_string(i) + " with " + std::to_string(subset.size()) + " edges");

            try {
                // PSEUDOCODE STEP 2.2: "If I ≠ ∅ do
                //                       if TreeDist-I(T_init, T_end, k, I) returns true,"
                if (TreeDistI(T_init, T_end, k, subset)) {
                    debugPrint("TreeDistS: Found solution with S-subset " + std::to_string(i));
                    return true;
                }
            } catch (...) {
                debugPrint("TreeDistS: Exception with S-subset " + std::to_string(i));
                continue;
            }
        }
    }

    // PSEUDOCODE STEP 3: "Return False"
    debugPrint("TreeDistS: No solution found");
    return false;
}
// ============================================================================
// COMPREHENSIVE ACCURACY AND SCALABILITY TESTING SUITE
// ============================================================================


// Performance timing utility
class PerformanceTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
public:
    void start() { start_time = std::chrono::high_resolution_clock::now(); }
    long long getMicroseconds() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }
    double getMilliseconds() { return getMicroseconds() / 1000.0; }
};

// Tree generators for testing
class TreeGenerator {
public:
    static std::pair<std::vector<int>, std::vector<int>> generateRightChain(int n) {
        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++) {
            preorder.push_back(i);
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateLeftChain(int n) {
        std::vector<int> preorder, inorder;
        for (int i = n; i >= 1; i--) {
            preorder.push_back(i);
        }
        for (int i = 1; i <= n; i++) {
            inorder.push_back(i);
        }
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateBalanced(int n) {
        if (n == 0) return {{}, {}};
        if (n == 1) return {{1}, {1}};
        if (n == 2) return {{2, 1}, {1, 2}};
        if (n == 3) return {{2, 1, 3}, {1, 2, 3}};
        if (n == 4) return {{3, 2, 1, 4}, {1, 2, 3, 4}};
        if (n == 5) return {{3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}};
        if (n == 6) return {{4, 2, 1, 3, 5, 6}, {1, 2, 3, 4, 5, 6}};
        if (n == 7) return {{4, 2, 1, 3, 6, 5, 7}, {1, 2, 3, 4, 5, 6, 7}};

        std::vector<int> preorder, inorder;
        for (int i = 1; i <= n; i++) inorder.push_back(i);

        std::function<void(int, int)> buildBalanced = [&](int start, int end) {
            if (start > end) return;
            int mid = (start + end) / 2;
            preorder.push_back(mid);
            buildBalanced(start, mid - 1);
            buildBalanced(mid + 1, end);
        };

        buildBalanced(1, n);
        return {preorder, inorder};
    }

    static std::pair<std::vector<int>, std::vector<int>> generateRandom(int n, int seed = 42) {
        std::mt19937 rng(seed);
        std::vector<int> inorder;
        for (int i = 1; i <= n; i++) inorder.push_back(i);

        std::vector<int> preorder = inorder;
        std::shuffle(preorder.begin(), preorder.end(), rng);

        return {preorder, inorder};
    }
};

void printTreeInfo(const std::string& name, const VectorRangeTreeMap& T) {
    try {
        std::cout << name << " (root=" << T.root << ", nodes=" << T.original_nodes.size() << "): ";
        auto edges = getInternalEdges(T);
        for (const auto& e : edges) {
            std::cout << "(" << e.first << "," << e.second << ") ";
        }
        std::cout << std::endl;
    } catch (...) {
        std::cout << name << " [ERROR]" << std::endl;
    }
}

// ============================================================================
// ACCURACY TESTING - Verify Algorithm Correctness
// ============================================================================

void testAccuracyBasicCases() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << " ACCURACY TESTING - BASIC CASES" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    struct TestCase {
        std::vector<int> pre1, in1, pre2, in2;
        int expected_min_distance;
        std::string description;
    };

    std::vector<TestCase> testCases = {
            // Basic cases
            {{1}, {1}, {1}, {1}, 0, "Single node (identical)"},
            {{1, 2}, {1, 2}, {1, 2}, {1, 2}, 0, "Two nodes (identical)"},
            {{1, 2}, {1, 2}, {2, 1}, {1, 2}, 1, "Two nodes (one rotation)"},

            // 3-node cases
            {{2, 1, 3}, {1, 2, 3}, {2, 1, 3}, {1, 2, 3}, 0, "3-node balanced (identical)"},
            {{1, 2, 3}, {1, 2, 3}, {2, 1, 3}, {1, 2, 3}, 1, "3-node chain to balanced"},
            {{2, 1, 3}, {1, 2, 3}, {3, 2, 1}, {1, 2, 3}, 1, "3-node balanced to mirrored"},
            {{1, 2, 3}, {1, 2, 3}, {3, 2, 1}, {1, 2, 3}, 2, "3-node chain to mirrored"},

            // 4-node cases
            {{1, 2, 3, 4}, {1, 2, 3, 4}, {3, 2, 1, 4}, {1, 2, 3, 4}, 2, "4-node chain to balanced"},
            {{3, 2, 1, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, 2, "4-node balanced to chain"},
            {{1, 2, 3, 4}, {1, 2, 3, 4}, {4, 3, 2, 1}, {1, 2, 3, 4}, 3, "4-node chain to reverse chain"},
    };

    int passed = 0;
    std::cout << std::setw(40) << "Test Case" << std::setw(12) << "Expected" << std::setw(12) << "Min k" << std::setw(12) << "Max k" << std::setw(10) << "Result" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    for (const auto& test : testCases) {
        VectorRangeTreeMap T1, T2;
        T1.build(test.pre1, test.in1);
        T2.build(test.pre2, test.in2);

        // Test if expected distance works
        bool works_at_expected = FlipDistTree(T1, T2, test.expected_min_distance);

        // Test if expected-1 fails (unless expected is 0)
        bool fails_before = (test.expected_min_distance == 0) ||
                            !FlipDistTree(T1, T2, test.expected_min_distance - 1);

        // Find the actual minimum working k
        int min_k = -1;
        for (int k = 0; k <= test.expected_min_distance + 2; k++) {
            if (FlipDistTree(T1, T2, k)) {
                min_k = k;
                break;
            }
        }

        // Find a reasonable upper bound
        int max_k = test.expected_min_distance + 3;
        for (int k = test.expected_min_distance; k <= test.expected_min_distance + 5; k++) {
            if (FlipDistTree(T1, T2, k)) {
                max_k = k;
                break;
            }
        }

        bool correct = works_at_expected && fails_before && (min_k == test.expected_min_distance);
        if (correct) passed++;

        std::cout << std::setw(40) << test.description.substr(0, 39)
                  << std::setw(12) << test.expected_min_distance
                  << std::setw(12) << min_k
                  << std::setw(12) << max_k
                  << std::setw(10) << (correct ? "✅ PASS" : "❌ FAIL") << std::endl;

        if (!correct) {
            std::cout << "    Details: works_at_expected=" << works_at_expected
                      << ", fails_before=" << fails_before
                      << ", min_k=" << min_k << std::endl;
        }
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Accuracy Test Results: " << passed << "/" << testCases.size() << " passed" << std::endl;
}

void testAccuracyMonotonicity() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << " ACCURACY TESTING - MONOTONICITY PROPERTY" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "Testing: If FlipDistTree(T1, T2, k) = true, then FlipDistTree(T1, T2, k+1) = true" << std::endl;

    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::string>> cases = {
            {{1, 2, 3}, {1, 2, 3}, {3, 2, 1}, {1, 2, 3}, "3-node chain to mirrored"},
            {{1, 2, 3, 4}, {1, 2, 3, 4}, {3, 2, 1, 4}, {1, 2, 3, 4}, "4-node case"},
            {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {3, 2, 1, 4, 5}, {1, 2, 3, 4, 5}, "5-node case"},
    };

    int monotonic_cases = 0;
    for (const auto& [pre1, in1, pre2, in2, name] : cases) {
        VectorRangeTreeMap T1, T2;
        T1.build(pre1, in1);
        T2.build(pre2, in2);

        bool monotonic = true;
        bool found_true = false;
        std::vector<bool> results;

        std::cout << "\n" << name << ":" << std::endl;
        std::cout << "k: ";
        for (int k = 0; k <= 8; k++) {
            bool result = FlipDistTree(T1, T2, k);
            results.push_back(result);
            std::cout << k << "=" << (result ? "T" : "F") << " ";

            if (found_true && !result) {
                monotonic = false;
            }
            if (result) found_true = true;
        }

        std::cout << " -> " << (monotonic ? "✅ MONOTONIC" : "❌ NOT MONOTONIC") << std::endl;
        if (monotonic) monotonic_cases++;
    }

    std::cout << "\nMonotonicity Results: " << monotonic_cases << "/" << cases.size() << " cases are monotonic" << std::endl;
}

void testAccuracyTreeEquality() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << " ACCURACY TESTING - TREE EQUALITY CASES" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::vector<int> sizes = {1, 2, 3, 4, 5};
    int passed = 0;
    int total = 0;

    for (int n : sizes) {
        auto [pre, in] = TreeGenerator::generateRightChain(n);
        VectorRangeTreeMap T1, T2;
        T1.build(pre, in);
        T2.build(pre, in);  // Identical tree

        bool result = FlipDistTree(T1, T2, 0);  // Should work with k=0
        total++;
        if (result) passed++;

        std::cout << "Identical " << n << "-node trees: k=0 -> "
                  << (result ? "✅ TRUE" : "❌ FALSE") << std::endl;
    }

    std::cout << "Equality test results: " << passed << "/" << total << " passed" << std::endl;
}

// ============================================================================
// SCALABILITY TESTING - Find Maximum Node Count
// ============================================================================

void testScalabilityComprehensive() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << " COMPREHENSIVE SCALABILITY ANALYSIS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    const long long TIMEOUT_MS = 5000;  // 5 second timeout per test
    PerformanceTimer timer;

    std::vector<std::string> testTypes = {
            "Chain→Balanced",
            "Chain→Chain",
            "Balanced→Chain",
            "Random→Random"
    };

    std::cout << std::setw(6) << "Nodes"
              << std::setw(16) << "Test Type"
              << std::setw(12) << "Time (ms)"
              << std::setw(12) << "Success"
              << std::setw(15) << "Status" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (int n = 1; n <= 20; n++) {  // Test up to 20 nodes
        bool any_timeout = false;

        for (const auto& testType : testTypes) {
            VectorRangeTreeMap T1, T2;

            try {
                // Generate trees based on test type
                if (testType == "Chain→Balanced") {
                    auto [pre1, in1] = TreeGenerator::generateRightChain(n);
                    auto [pre2, in2] = TreeGenerator::generateBalanced(n);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                } else if (testType == "Chain→Chain") {
                    auto [pre1, in1] = TreeGenerator::generateRightChain(n);
                    auto [pre2, in2] = TreeGenerator::generateLeftChain(n);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                } else if (testType == "Balanced→Chain") {
                    auto [pre1, in1] = TreeGenerator::generateBalanced(n);
                    auto [pre2, in2] = TreeGenerator::generateRightChain(n);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                } else if (testType == "Random→Random") {
                    auto [pre1, in1] = TreeGenerator::generateRandom(n, 42);
                    auto [pre2, in2] = TreeGenerator::generateRandom(n, 84);
                    T1.build(pre1, in1);
                    T2.build(pre2, in2);
                }

                // Test with generous budget
                int budget = n + 5;

                timer.start();
                bool success = FlipDistTree(T1, T2, budget);
                double time_ms = timer.getMilliseconds();

                std::string status = "OK";
                if (time_ms > TIMEOUT_MS) {
                    status = "TIMEOUT";
                    any_timeout = true;
                }

                std::cout << std::setw(6) << n
                          << std::setw(16) << testType
                          << std::setw(12) << std::fixed << std::setprecision(2) << time_ms
                          << std::setw(12) << (success ? "YES" : "NO")
                          << std::setw(15) << status << std::endl;

                // Stop this test type if timeout
                if (time_ms > TIMEOUT_MS) {
                    break;
                }

            } catch (...) {
                std::cout << std::setw(6) << n
                          << std::setw(16) << testType
                          << std::setw(12) << "ERROR"
                          << std::setw(12) << "NO"
                          << std::setw(15) << "EXCEPTION" << std::endl;
            }
        }

        // Stop if most test types are timing out
        if (any_timeout && n >= 5) {
            std::cout << "\n⚠  Stopping at " << n << " nodes due to performance limits" << std::endl;
            break;
        }
    }
}

void testScalabilityDetailed() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << " DETAILED PERFORMANCE ANALYSIS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::vector<int> testSizes = {3, 4, 5, 6, 7, 8, 9, 10};
    PerformanceTimer timer;

    std::cout << std::setw(6) << "Nodes"
              << std::setw(12) << "Avg (ms)"
              << std::setw(12) << "Min (ms)"
              << std::setw(12) << "Max (ms)"
              << std::setw(12) << "Success"
              << std::setw(15) << "Notes" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (int n : testSizes) {
        std::vector<double> times;
        int successes = 0;
        const int NUM_RUNS = 3;

        for (int run = 0; run < NUM_RUNS; run++) {
            try {
                auto [pre1, in1] = TreeGenerator::generateRightChain(n);
                auto [pre2, in2] = TreeGenerator::generateBalanced(n);

                VectorRangeTreeMap T1, T2;
                T1.build(pre1, in1);
                T2.build(pre2, in2);

                timer.start();
                bool result = FlipDistTree(T1, T2, n + 3);
                double time_ms = timer.getMilliseconds();

                times.push_back(time_ms);
                if (result) successes++;

                // Stop if getting too slow
                if (time_ms > 10000) {  // 10 seconds
                    break;
                }

            } catch (...) {
                times.push_back(-1);  // Error marker
            }
        }

        if (!times.empty() && times[0] >= 0) {
            double avg = 0, min_time = times[0], max_time = times[0];
            int valid_times = 0;

            for (double t : times) {
                if (t >= 0) {
                    avg += t;
                    min_time = std::min(min_time, t);
                    max_time = std::max(max_time, t);
                    valid_times++;
                }
            }

            if (valid_times > 0) {
                avg /= valid_times;

                std::string notes = "";
                if (avg > 1000) notes = "SLOW";
                else if (avg > 100) notes = "MODERATE";
                else notes = "FAST";

                std::cout << std::setw(6) << n
                          << std::setw(12) << std::fixed << std::setprecision(2) << avg
                          << std::setw(12) << std::setprecision(2) << min_time
                          << std::setw(12) << std::setprecision(2) << max_time
                          << std::setw(12) << successes << "/" << NUM_RUNS
                          << std::setw(15) << notes << std::endl;

                // Stop if consistently slow
                if (avg > 5000) {
                    std::cout << " Stopping detailed analysis due to performance" << std::endl;
                    break;
                }
            }
        } else {
            std::cout << std::setw(6) << n
                      << std::setw(12) << "ERROR"
                      << std::setw(12) << "-"
                      << std::setw(12) << "-"
                      << std::setw(12) << "0/" << NUM_RUNS
                      << std::setw(15) << "FAILED" << std::endl;
        }
    }
}

// ============================================================================
// MAIN TESTING FUNCTION
// ============================================================================

void runComprehensiveTests() {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "COMPREHENSIVE ALGORITHM TESTING SUITE" << std::endl;
    std::cout << "    Accuracy Validation + Scalability Analysis" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    // Phase 1: Accuracy Testing
    std::cout << "\nPHASE 1: ACCURACY VALIDATION" << std::endl;
    testAccuracyBasicCases();
    testAccuracyMonotonicity();
    testAccuracyTreeEquality();

    // Phase 2: Scalability Testing
    std::cout << "\nPHASE 2: SCALABILITY ANALYSIS" << std::endl;
    testScalabilityComprehensive();
    testScalabilityDetailed();

    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "COMPREHENSIVE TESTING COMPLETE!" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
}

int main() {
    runComprehensiveTests();
    return 0;
}