// Conflict detection and edge utilities for FlipDist
#pragma once

#include <vector>
#include <utility>
#include "../rotation_tree.h"
#include "types.h"

// All internal parent->child edges in the tree.
std::vector<std::pair<int,int>> getInternalEdges(const VectorRangeTreeMap &T);

// Incident edges (parent/children) for a node.
std::vector<std::pair<int,int>> getIncidentEdges(const VectorRangeTreeMap &T, int node);

// Count internal edges.
int countInternalEdges(const VectorRangeTreeMap &T);

// Lower bound on rotations based on edge differences.
int lowerBoundEdgeDifference(const VectorRangeTreeMap &A,
                             const VectorRangeTreeMap &B);

// Build the list of conflicting diagonals (present in start but not target).
std::vector<DiagonalEdge> collectConflictingEdges(
        const VectorRangeTreeMap &start,
        const VectorRangeTreeMap &target);

// Greedy maximal independent set of conflicting edges.
std::vector<std::pair<int,int>> buildMaxIndependentSet(
        const VectorRangeTreeMap &start,
        const std::vector<DiagonalEdge> &conflicts);
