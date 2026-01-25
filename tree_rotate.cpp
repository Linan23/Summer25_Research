// Houses the left/right rotation implementations used during BFS expansions
// and demo path reconstruction. Rotations reuse the metadata maintenance
// hooks defined in tree_core.cpp and keep the tree’s range and parent tables
// consistent after each adjustment.

#include "rotation_tree.h"

// Performs a standard AVL-style left rotation at node x and fixes metadata.
void VectorRangeTreeMap::rotateLeft(int x) {
    int y = getRightChild(x);
    if (y == NO_CHILD) return;

    int yl = getLeftChild(y);
    int xp = getParent(x);

    // rewire
    parents[y] = xp;
    if (xp != NO_PARENT) {
        if (edges[xp].first == x)  edges[xp].first  = y;
        else                       edges[xp].second = y;
    } else {
        root = y;
    }
    parents[x] = y;
    edges[y].first = x;

    edges[x].second = yl;
    if (yl >= 0) parents[yl] = x;

    // update ranges bottom-up
    updateNodeRange(x);
    updateNodeRange(y);
    recomputeUpwardsFrom(getParent(y));
    invalidateSignature();

#ifndef NDEBUG
    verify();
#endif
}

// Performs a standard AVL-style right rotation at node x and fixes metadata.
void VectorRangeTreeMap::rotateRight(int x) {
    int y = getLeftChild(x);
    if (y == NO_CHILD) return;

    int yr = getRightChild(y);
    int xp = getParent(x);

    // rewire
    parents[y] = xp;
    if (xp != NO_PARENT) {
        if (edges[xp].first == x)  edges[xp].first  = y;
        else                       edges[xp].second = y;
    } else {
        root = y;
    }
    parents[x] = y;
    edges[y].second = x;

    edges[x].first = yr;
    if (yr >= 0) parents[yr] = x;

    // update ranges bottom-up
    updateNodeRange(x);
    updateNodeRange(y);
    recomputeUpwardsFrom(getParent(y));
    invalidateSignature();

#ifndef NDEBUG
    verify();
#endif
}
