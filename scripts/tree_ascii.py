#!/usr/bin/env python3
"""
Render binary trees from the canonical traversal format used in the project.

Usage examples:

  # Visualise a simple random tree (seeded for determinism)
  python3 scripts/tree_ascii.py --random 6 --seed 42

  # Visualise the left and right combs that form the hard cases
  python3 scripts/tree_ascii.py --comb 6 --lean left
  python3 scripts/tree_ascii.py --comb 6 --lean right

  # Visualise a tree directly from the canonical string stored in the CSV logs
  python3 scripts/tree_ascii.py --tree "P:12,11,10,9,8,7,6,5,4,3,2,1;I:1,2,3,4,5,6,7,8,9,10,11,12"
"""

import argparse
import random
import sys
from typing import Dict, List, Optional, Tuple


class Node:
    __slots__ = ("value", "left", "right")

    def __init__(self, value: int) -> None:
        self.value = value
        self.left: Optional["Node"] = None
        self.right: Optional["Node"] = None


def parse_canonical(tree: str) -> Tuple[List[int], List[int]]:
    tree = tree.strip()
    if not tree.startswith("P:") or ";I:" not in tree:
        raise ValueError("Expected canonical string format 'P:...;I:...'")
    pre_part, in_part = tree[2:].split(";I:", 1)
    preorder = [int(tok) for tok in pre_part.split(",") if tok]
    inorder = [int(tok) for tok in in_part.split(",") if tok]
    if len(preorder) != len(inorder):
        raise ValueError("Preorder and inorder lengths differ")
    return preorder, inorder


def build_from_traversals(
    preorder: List[int], inorder: List[int]
) -> Optional[Node]:
    index: Dict[int, int] = {value: pos for pos, value in enumerate(inorder)}

    def helper(ps: int, pe: int, is_: int, ie: int) -> Optional[Node]:
        if ps >= pe or is_ >= ie:
            return None
        root_val = preorder[ps]
        root_pos = index[root_val]
        node = Node(root_val)
        left_len = root_pos - is_
        node.left = helper(ps + 1, ps + 1 + left_len, is_, root_pos)
        node.right = helper(ps + 1 + left_len, pe, root_pos + 1, ie)
        return node

    return helper(0, len(preorder), 0, len(inorder))


def make_comb(n: int, lean_right: bool) -> Tuple[List[int], List[int]]:
    if n <= 0:
        return [], []

    inorder = list(range(1, n + 1))
    if lean_right:
        preorder = list(range(1, n + 1))
    else:
        preorder = list(range(n, 0, -1))
    return preorder, inorder


def make_random(n: int, seed: int) -> Tuple[List[int], List[int]]:
    rng = random.Random(seed)

    def build(size: int) -> Optional[Node]:
        if size <= 0:
            return None
        root = Node(-1)
        left_size = rng.randint(0, size - 1)
        right_size = size - 1 - left_size
        root.left = build(left_size)
        root.right = build(right_size)
        return root

    root = build(n)
    preorder: List[int] = []
    inorder: List[int] = []

    def assign_inorder(node: Optional[Node], next_id: List[int]) -> None:
        if not node:
            return
        assign_inorder(node.left, next_id)
        node.value = next_id[0]
        next_id[0] += 1
        inorder.append(node.value)
        assign_inorder(node.right, next_id)

    def record_preorder(node: Optional[Node]) -> None:
        if not node:
            return
        preorder.append(node.value)
        record_preorder(node.left)
        record_preorder(node.right)

    assign_inorder(root, [1])
    record_preorder(root)
    return preorder, inorder


def render_ascii(root: Optional[Node]) -> List[str]:
    lines: List[str] = []

    def walk(node: Optional[Node], prefix: str, is_left: bool) -> None:
        if not node:
            return
        connector = "└── " if is_left else "┌── "
        walk(node.right, prefix + ("│   " if is_left else "    "), False)
        lines.append(prefix + connector + str(node.value))
        walk(node.left, prefix + ("    " if is_left else "│   "), True)

    if root:
        walk(root, "", True)
    else:
        lines.append("(empty)")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description="Render binary trees as ASCII art")
    parser.add_argument("--tree", type=str, help="Canonical tree string (P:...;I:...)")
    parser.add_argument("--comb", type=int, help="Generate comb tree with N internal nodes")
    parser.add_argument("--lean", choices=["left", "right"], default="right",
                        help="Lean direction for comb generation")
    parser.add_argument("--random", type=int, help="Generate random tree with N nodes")
    parser.add_argument("--seed", type=int, default=2025, help="Seed for random generator")
    args = parser.parse_args()

    if args.tree:
        preorder, inorder = parse_canonical(args.tree)
    elif args.comb is not None:
        preorder, inorder = make_comb(args.comb, lean_right=(args.lean == "right"))
    elif args.random is not None:
        preorder, inorder = make_random(args.random, args.seed)
    else:
        parser.error("Provide --tree, --comb, or --random")
        return

    root = build_from_traversals(preorder, inorder)
    canonical = f"P:{','.join(map(str, preorder))};I:{','.join(map(str, inorder))}"

    print("Canonical:", canonical)
    print("\nASCII Visualisation:\n")
    for line in render_ascii(root):
        print(line)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover - quick CLI guard
        sys.stderr.write(f"Error: {exc}\n")
        sys.exit(1)
