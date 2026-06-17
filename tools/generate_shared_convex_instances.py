#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
import random


@dataclass
class Node:
    value: int
    left: "Node | None" = None
    right: "Node | None" = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate shared convex-polygon benchmark instances for FlipDist and A*"
    )
    parser.add_argument(
        "--case",
        choices=["random", "simple", "comb"],
        default="random",
        help="'simple' is the clearer alias for the older 'comb'.",
    )
    parser.add_argument("--n", type=int, required=True, help="Number of BST nodes / number of triangles")
    parser.add_argument("--seed-min", type=int, default=0, help="First seed to generate")
    parser.add_argument("--seed-max", type=int, default=0, help="Last seed to generate")
    parser.add_argument(
        "--output-dir",
        "--output",
        dest="output_dir",
        required=True,
        help="Root directory for generated instances",
    )
    parser.add_argument(
        "--flipdist-binary",
        default="./build/flipdist",
        help="Path recorded in the manifest for FlipDist runs",
    )
    parser.add_argument(
        "--astar-binary",
        default="third_party/AStarFlipDistance/build-nogurobi/A_star_for_flipdistance",
        help="Path recorded in the manifest for A* runs",
    )
    return parser.parse_args()


def insert_bst(root: Node | None, value: int) -> Node:
    if root is None:
        return Node(value)
    cur = root
    while True:
        if value < cur.value:
            if cur.left is None:
                cur.left = Node(value)
                break
            cur = cur.left
        else:
            if cur.right is None:
                cur.right = Node(value)
                break
            cur = cur.right
    return root


def preorder(node: Node | None, out: list[int]) -> None:
    if node is None:
        return
    out.append(node.value)
    preorder(node.left, out)
    preorder(node.right, out)


def inorder(node: Node | None, out: list[int]) -> None:
    if node is None:
        return
    inorder(node.left, out)
    out.append(node.value)
    inorder(node.right, out)


def generate_random_tree(n: int, seed: int) -> tuple[list[int], list[int]]:
    values = list(range(1, n + 1))
    rng = random.Random(seed)
    rng.shuffle(values)
    root: Node | None = None
    for value in values:
        root = insert_bst(root, value)
    pre: list[int] = []
    ino: list[int] = []
    preorder(root, pre)
    inorder(root, ino)
    return pre, ino


def generate_left_chain(n: int) -> tuple[list[int], list[int]]:
    return list(range(n, 0, -1)), list(range(1, n + 1))


def generate_right_chain(n: int) -> tuple[list[int], list[int]]:
    seq = list(range(1, n + 1))
    return seq[:], seq


def build_from_traversals(pre: list[int], ino: list[int]) -> Node | None:
    if not pre:
        return None
    index = {value: pos for pos, value in enumerate(ino)}

    def rec(ps: int, pe: int, is_: int, ie: int) -> Node | None:
        if ps >= pe or is_ >= ie:
            return None
        root_value = pre[ps]
        pos = index[root_value]
        left_size = pos - is_
        node = Node(root_value)
        node.left = rec(ps + 1, ps + 1 + left_size, is_, pos)
        node.right = rec(ps + 1 + left_size, pe, pos + 1, ie)
        return node

    return rec(0, len(pre), 0, len(ino))


def tree_to_canonical(pre: list[int], ino: list[int]) -> str:
    return "P:" + ",".join(map(str, pre)) + ";I:" + ",".join(map(str, ino))


def triangles_from_tree(pre: list[int], ino: list[int]) -> list[tuple[int, int, int]]:
    root = build_from_traversals(pre, ino)
    if root is None:
        return []
    n = len(ino)
    triangles: list[tuple[int, int, int]] = []

    def rec(node: Node | None, left_boundary: int, right_boundary: int) -> None:
        if node is None:
            return
        tri = tuple(sorted((left_boundary, node.value, right_boundary)))
        triangles.append(tri)
        rec(node.left, left_boundary, node.value)
        rec(node.right, node.value, right_boundary)

    rec(root, 0, n + 1)
    return sorted(triangles)


def convex_points(vertex_count: int) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    for i in range(vertex_count):
        theta = (2.0 * math.pi * i) / vertex_count
        x = 1000.0 * math.cos(theta)
        y = 700.0 * math.sin(theta)
        points.append((x, y))
    return points


def validate_triangulation(n: int, triangles: list[tuple[int, int, int]]) -> None:
    expected_triangles = n
    if len(triangles) != expected_triangles:
        raise ValueError(f"Expected {expected_triangles} triangles, found {len(triangles)}")

    vertex_count = n + 2
    boundary_edges = {
        tuple(sorted((i, i + 1))) for i in range(vertex_count - 1)
    }
    boundary_edges.add((0, vertex_count - 1))

    edge_counts: dict[tuple[int, int], int] = {}
    for a, b, c in triangles:
        for edge in (tuple(sorted((a, b))), tuple(sorted((b, c))), tuple(sorted((a, c)))):
            edge_counts[edge] = edge_counts.get(edge, 0) + 1

    for edge, count in edge_counts.items():
        if edge in boundary_edges:
            if count != 1:
                raise ValueError(f"Boundary edge {edge} appears {count} times")
        else:
            if count != 2:
                raise ValueError(f"Internal edge {edge} appears {count} times")

    expected_internal = n - 1
    actual_internal = sum(1 for edge in edge_counts if edge not in boundary_edges)
    if actual_internal != expected_internal:
        raise ValueError(f"Expected {expected_internal} internal edges, found {actual_internal}")


def write_points(path: Path, points: list[tuple[float, float]]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for x, y in points:
            fh.write(f"{x:.6f} {y:.6f}\n")


def write_triangles(path: Path, triangles: list[tuple[int, int, int]]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for a, b, c in triangles:
            fh.write(f"{a} {b} {c}\n")


def generate_case(case_type: str, n: int, seed: int) -> tuple[tuple[list[int], list[int]], tuple[list[int], list[int]]]:
    if case_type == "random":
        return generate_random_tree(n, seed * 2), generate_random_tree(n, seed * 2 + 1)
    return generate_left_chain(n), generate_right_chain(n)


def main() -> int:
    cfg = parse_args()
    if cfg.case == "comb":
        cfg.case = "simple"
    if cfg.n <= 0:
        raise SystemExit("--n must be positive")
    if cfg.seed_max < cfg.seed_min:
        raise SystemExit("--seed-max must be >= --seed-min")

    root = Path(cfg.output_dir).resolve()
    vertex_count = cfg.n + 2
    points = convex_points(vertex_count)

    generated = 0
    for seed in range(cfg.seed_min, cfg.seed_max + 1):
        if cfg.seed_min == cfg.seed_max:
            case_dir = root
        else:
            case_dir = root / f"n_{cfg.n}" / f"set_{seed}"
        tri_dir = case_dir / "triangulation"
        tri_dir.mkdir(parents=True, exist_ok=True)

        (pre_a, in_a), (pre_b, in_b) = generate_case(cfg.case, cfg.n, seed)
        tree_a = tree_to_canonical(pre_a, in_a)
        tree_b = tree_to_canonical(pre_b, in_b)
        triangles_a = triangles_from_tree(pre_a, in_a)
        triangles_b = triangles_from_tree(pre_b, in_b)
        validate_triangulation(cfg.n, triangles_a)
        validate_triangulation(cfg.n, triangles_b)

        points_path = case_dir / "points"
        tri_a_path = tri_dir / "t_0"
        tri_b_path = tri_dir / "t_1"
        tree_a_path = case_dir / "tree_a.txt"
        tree_b_path = case_dir / "tree_b.txt"
        manifest_path = case_dir / "manifest.json"

        write_points(points_path, points)
        write_triangles(tri_a_path, triangles_a)
        write_triangles(tri_b_path, triangles_b)
        tree_a_path.write_text(tree_a + "\n", encoding="utf-8")
        tree_b_path.write_text(tree_b + "\n", encoding="utf-8")

        manifest = {
            "case_type": cfg.case,
            "n": cfg.n,
            "seed": seed,
            "vertex_count": vertex_count,
            "tree_a": tree_a,
            "tree_b": tree_b,
            "points_path": str(points_path),
            "triangulation_1_path": str(tri_a_path),
            "triangulation_2_path": str(tri_b_path),
            "flipdist_tree_a_file": str(tree_a_path),
            "flipdist_tree_b_file": str(tree_b_path),
            "flipdist_command": [
                cfg.flipdist_binary,
                "--tree-a-file",
                str(tree_a_path),
                "--tree-b-file",
                str(tree_b_path),
                "--max-k",
                str(max(1, 3 * cfg.n + 10)),
            ],
            "astar_command": [
                cfg.astar_binary,
                str(points_path),
                str(tri_a_path),
                str(tri_b_path),
                "simple",
                str(case_dir),
                "astar_simple.json",
            ],
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        generated += 1

        print(f"Generated shared convex instance: n={cfg.n} seed={seed} -> {case_dir}")

    print(f"Generated {generated} shared convex instance(s) under {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
