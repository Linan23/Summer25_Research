// Basic shared types for FlipDist
#pragma once

#include <utility>

struct DiagonalEdge
{
    std::pair<int,int> diag; // polygon endpoints (L,R), L < R
    int parent;              // oriented parent node in the current tree
    int child;               // oriented child node in the current tree
};
