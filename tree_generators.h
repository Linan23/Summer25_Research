#ifndef TREE_GENERATORS_H
#define TREE_GENERATORS_H

#include <random>
#include <vector>

struct Traversals {
    std::vector<int> preorder;
    std::vector<int> inorder;
};

Traversals makeCombTraversals(int n, bool rightComb);
Traversals makeRandomTraversals(int n, std::mt19937& rng);

#endif // TREE_GENERATORS_H
