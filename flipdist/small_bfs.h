// Fast exact rotation distance lookup for very small trees (used to accelerate local BFS fallbacks).
#pragma once

#include "../rotation_tree.h"

namespace flipdist {

// Returns the exact rotation distance between A and B when both trees have <= 9 nodes,
// or -1 if the distance exceeds cap. Returns -2 when the lookup is not applicable.
int smallRotationDistance(const VectorRangeTreeMap &A, const VectorRangeTreeMap &B, int cap);

} // namespace flipdist
