/*
  Short terminology used by the browser visualizer.
  These definitions are intentionally simple so a new reader can inspect a
  solver run before reading the research notes.
*/
export const TERMS = [
  {
    term: "Tree A / Tree B",
    text: "The two rooted binary trees in one solver case."
  },
  {
    term: "Rotation",
    text: "A local tree move that changes one parent-child relationship while keeping the inorder node order."
  },
  {
    term: "Flip distance",
    text: "The minimum number of rotations needed to transform one tree into the other."
  },
  {
    term: "Directed result",
    text: "The solver reports both a->b and b->a. The distance should agree, but runtime can differ."
  },
  {
    term: "Exact distance",
    text: "The distance reported by the C++ solver when the row has status ok."
  },
  {
    term: "Simple case",
    text: "A case with immediate structural progress, such as a common edge or a directly helpful move."
  },
  {
    term: "Complex case",
    text: "A case without immediate simple progress, so the recursive search usually does more work."
  },
  {
    term: "Visual path",
    text: "A browser-reconstructed shortest rotation sequence. If the exact shortest path is too large to reconstruct, animation is disabled."
  }
];
