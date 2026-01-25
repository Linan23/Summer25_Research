// Miscellaneous helpers shared across FlipDist implementation.
#pragma once

#include <string>
#include <vector>
#include <utility>
#include <unordered_set>
#include "types.h"

struct PairHash;
struct PairEq;

std::string canonicalEdgeListKey(const std::vector<std::pair<int,int>> &edges);
std::string canonicalDiagonalPairListKey(const std::vector<std::pair<DiagonalEdge, DiagonalEdge>> &pairs);
std::string canonicalDiagonalSetKey(const std::unordered_set<std::pair<int,int>, PairHash, PairEq> &diags);
std::string makeDiagonalPairKey(const DiagonalEdge &a, const DiagonalEdge &b);
std::string makeRangeMemoKey(const std::string &treePairSig,
                             const std::string &conflictSig,
                             const std::pair<int,int> &range);
