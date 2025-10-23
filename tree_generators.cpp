// Helper factories for producing random or comb-shaped traversal pairs. The
// demo harness, probes, and comparison script all rely on these helpers so we
// consistently build identical `VectorRangeTreeMap` inputs on both the C++ and
// Java sides.

#include "tree_generators.h"

#include <algorithm>
#include <numeric>

namespace {

struct Shape {
    Shape* L{nullptr};
    Shape* R{nullptr};
    int label{-1};
};

Shape* buildRandomShape(int n, std::mt19937& rng) {
    if (n <= 0) return nullptr;
    if (n == 1) return new Shape();

    std::uniform_int_distribution<int> pick(0, n - 1);
    int leftSize = pick(rng);
    int rightSize = n - 1 - leftSize;

    Shape* node = new Shape();
    node->L = buildRandomShape(leftSize, rng);
    node->R = buildRandomShape(rightSize, rng);
    return node;
}

void assignInorder(Shape* node, int& next) {
    if (!node) return;
    assignInorder(node->L, next);
    node->label = next++;
    assignInorder(node->R, next);
}

void dumpInorder(Shape* node, std::vector<int>& out) {
    if (!node) return;
    dumpInorder(node->L, out);
    out.push_back(node->label);
    dumpInorder(node->R, out);
}

void dumpPreorder(Shape* node, std::vector<int>& out) {
    if (!node) return;
    out.push_back(node->label);
    dumpPreorder(node->L, out);
    dumpPreorder(node->R, out);
}

void freeShape(Shape* node) {
    if (!node) return;
    freeShape(node->L);
    freeShape(node->R);
    delete node;
}

} // namespace

Traversals makeCombTraversals(int n, bool rightComb) {
    Traversals t;
    if (n <= 0) return t;

    t.inorder.resize(n);
    std::iota(t.inorder.begin(), t.inorder.end(), 1);

    t.preorder.resize(n);
    if (rightComb) {
        std::iota(t.preorder.begin(), t.preorder.end(), 1);
    } else {
        for (int i = 0; i < n; ++i) t.preorder[i] = n - i;
    }
    return t;
}

Traversals makeRandomTraversals(int n, std::mt19937& rng) {
    Traversals t;
    if (n <= 0) return t;

    Shape* root = buildRandomShape(n, rng);
    int next = 1;
    assignInorder(root, next);
    dumpInorder(root, t.inorder);
    dumpPreorder(root, t.preorder);
    freeShape(root);
    return t;
}
