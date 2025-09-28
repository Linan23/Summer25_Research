#include "A_tree.h"

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

#ifndef NDEBUG
    verify();
#endif
}

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

#ifndef NDEBUG
    verify();
#endif
}
