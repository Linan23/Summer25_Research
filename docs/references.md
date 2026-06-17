# References

## Primary Algorithmic Reference

This repository builds on the exact fixed-parameter approach of:

Haohong Li and Ge Xia. "An O(3.82^k) Time FPT Algorithm for Convex Flip Distance." In *40th International Symposium on Theoretical Aspects of Computer Science (STACS 2023)*, Article 44. DOI: `10.4230/LIPIcs.STACS.2023.44`.

Links:

- Published PDF: https://d-nb.info/1367146364/34
- DOI landing page: https://doi.org/10.4230/LIPIcs.STACS.2023.44
- Full version: https://arxiv.org/abs/2209.13134

## Problem Context

The paper studies Convex Flip Distance: given two triangulations of a convex polygon and a parameter `k`, decide whether their flip distance is at most `k`. It also uses the standard equivalence between convex polygon triangulations and rooted full binary trees: an edge flip in a triangulation corresponds to a rotation in the associated binary tree.

The implementation in this repository works in the binary-tree representation while preserving the exact flip-distance semantics and the Li-Xia search structure used by the paper.
