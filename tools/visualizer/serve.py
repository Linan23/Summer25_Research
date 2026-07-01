#!/usr/bin/env python3
"""Serve the FlipDist browser visualizer from a local static server.

The visualizer has no frontend build step. This launcher also exposes a small
local API that runs the existing flipdist binary to generate cases on demand.
"""

from __future__ import annotations

import argparse
import functools
import http.server
import json
import socket
import socketserver
import subprocess
import sys
import tempfile
import urllib.parse
import webbrowser
from pathlib import Path


class ReusableThreadingServer(socketserver.ThreadingTCPServer):
    allow_reuse_address = True


class VisualizerHandler(http.server.SimpleHTTPRequestHandler):
    """Static-file handler with a local flipdist generation endpoint."""

    def __init__(self, *args, directory: str, repo_root: Path, solver: Path, solver_timeout: float, **kwargs):
        self.repo_root = repo_root
        self.solver = solver
        self.solver_timeout = solver_timeout
        super().__init__(*args, directory=directory, **kwargs)

    def do_GET(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path == "/api/generate":
            self.handle_generate(parsed.query)
            return
        super().do_GET()

    def send_json(self, status: int, payload: dict) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def parse_int(self, query: dict[str, list[str]], name: str, default: int,
                  min_value: int, max_value: int) -> int:
        raw = query.get(name, [str(default)])[0]
        try:
            value = int(raw)
        except ValueError as exc:
            raise ValueError(f"{name} must be an integer.") from exc
        if value < min_value or value > max_value:
            raise ValueError(f"{name} must be between {min_value} and {max_value}.")
        return value

    @staticmethod
    def parse_tree_encoding(encoding: str) -> tuple[list[int], list[int]]:
        parts: dict[str, str] = {}
        for raw in encoding.strip().split(";"):
            if ":" not in raw:
                continue
            key, value = raw.split(":", 1)
            parts[key.strip()] = value.strip()
        if "P" not in parts or "I" not in parts:
            raise ValueError("Tree encoding must contain P:...;I:...")

        def parse_list(text: str) -> list[int]:
            if not text:
                return []
            return [int(item.strip()) for item in text.split(",") if item.strip()]

        preorder = parse_list(parts["P"])
        inorder = parse_list(parts["I"])
        if len(preorder) != len(inorder):
            raise ValueError("Preorder and inorder lengths differ.")
        return preorder, inorder

    @staticmethod
    def encode_tree(preorder: list[int], inorder: list[int]) -> str:
        return "P:" + ",".join(str(x) for x in preorder) + ";I:" + ",".join(str(x) for x in inorder)

    @classmethod
    def build_tree(cls, encoding: str) -> dict:
        preorder, inorder = cls.parse_tree_encoding(encoding)
        pos = {value: idx for idx, value in enumerate(inorder)}
        nodes: dict[int, dict] = {}

        def build(ps: int, pe: int, is_: int, ie: int, parent: int | None) -> int | None:
            if ps > pe:
                return None
            value = preorder[ps]
            idx = pos[value]
            left_size = idx - is_
            nodes[value] = {"left": None, "right": None, "parent": parent}
            nodes[value]["left"] = build(ps + 1, ps + left_size, is_, idx - 1, value)
            nodes[value]["right"] = build(ps + left_size + 1, pe, idx + 1, ie, value)
            return value

        root = build(0, len(preorder) - 1, 0, len(inorder) - 1, None) if preorder else None
        return {"root": root, "nodes": nodes, "inorder": inorder}

    @staticmethod
    def clone_tree(tree: dict) -> dict:
        return {
            "root": tree["root"],
            "inorder": list(tree["inorder"]),
            "nodes": {
                value: {
                    "left": node["left"],
                    "right": node["right"],
                    "parent": node["parent"],
                }
                for value, node in tree["nodes"].items()
            },
        }

    @classmethod
    def tree_preorder(cls, tree: dict) -> list[int]:
        out: list[int] = []

        def visit(value: int | None) -> None:
            if value is None:
                return
            node = tree["nodes"].get(value)
            if node is None:
                return
            out.append(value)
            visit(node["left"])
            visit(node["right"])

        visit(tree["root"])
        return out

    @classmethod
    def tree_encoding(cls, tree: dict) -> str:
        return cls.encode_tree(cls.tree_preorder(tree), tree["inorder"])

    @staticmethod
    def rotate_left(tree: dict, pivot_value: int) -> bool:
        pivot = tree["nodes"].get(pivot_value)
        if not pivot or pivot["right"] is None:
            return False
        child_value = pivot["right"]
        child = tree["nodes"][child_value]
        beta_value = child["left"]
        parent_value = pivot["parent"]

        child["parent"] = parent_value
        if parent_value is None:
            tree["root"] = child_value
        else:
            parent = tree["nodes"][parent_value]
            if parent["left"] == pivot_value:
                parent["left"] = child_value
            else:
                parent["right"] = child_value

        child["left"] = pivot_value
        pivot["parent"] = child_value
        pivot["right"] = beta_value
        if beta_value is not None:
            tree["nodes"][beta_value]["parent"] = pivot_value
        return True

    @staticmethod
    def rotate_right(tree: dict, pivot_value: int) -> bool:
        pivot = tree["nodes"].get(pivot_value)
        if not pivot or pivot["left"] is None:
            return False
        child_value = pivot["left"]
        child = tree["nodes"][child_value]
        beta_value = child["right"]
        parent_value = pivot["parent"]

        child["parent"] = parent_value
        if parent_value is None:
            tree["root"] = child_value
        else:
            parent = tree["nodes"][parent_value]
            if parent["left"] == pivot_value:
                parent["left"] = child_value
            else:
                parent["right"] = child_value

        child["right"] = pivot_value
        pivot["parent"] = child_value
        pivot["left"] = beta_value
        if beta_value is not None:
            tree["nodes"][beta_value]["parent"] = pivot_value
        return True

    @classmethod
    def rotation_neighbors(cls, encoding: str) -> list[tuple[str, str]]:
        tree = cls.build_tree(encoding)
        out: list[tuple[str, str]] = []
        seen: set[str] = set()
        for value in tree["inorder"]:
            node = tree["nodes"].get(value)
            if node and node["right"] is not None:
                nxt = cls.clone_tree(tree)
                cls.rotate_left(nxt, value)
                enc = cls.tree_encoding(nxt)
                if enc not in seen:
                    seen.add(enc)
                    out.append((f"rotateLeft({value})", enc))
            if node and node["left"] is not None:
                nxt = cls.clone_tree(tree)
                cls.rotate_right(nxt, value)
                enc = cls.tree_encoding(nxt)
                if enc not in seen:
                    seen.add(enc)
                    out.append((f"rotateRight({value})", enc))
        return out

    def solve_custom_distance(self, source: str, target: str, max_k: int, tmpdir: Path) -> int | None:
        a_path = tmpdir / "a.tree"
        b_path = tmpdir / "b.tree"
        a_path.write_text(source, encoding="utf-8")
        b_path.write_text(target, encoding="utf-8")
        cmd = [
            str(self.solver),
            "--tree-a-file",
            str(a_path),
            "--tree-b-file",
            str(b_path),
            "--max-k",
            str(max_k),
            "--bfs-cap",
            "1",
        ]
        try:
            proc = subprocess.run(
                cmd,
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                check=False,
                timeout=self.solver_timeout,
            )
        except subprocess.TimeoutExpired:
            return None
        if proc.returncode != 0:
            return None
        for line in proc.stdout.splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("direction") == "a->b" and row.get("status") == "ok":
                return int(row.get("distance"))
        return None

    def build_witness(self, source: str, target: str, distance: int) -> dict:
        if distance < 0:
            return {"status": "unavailable", "message": "distance was not solved"}
        steps = [{"step": 0, "move": "start", "tree": source}]
        current = source
        remaining = distance

        with tempfile.TemporaryDirectory(prefix="flipdist_witness_") as tmp:
            tmpdir = Path(tmp)
            while remaining > 0:
                found: tuple[str, str] | None = None
                for move, candidate in self.rotation_neighbors(current):
                    dist = self.solve_custom_distance(candidate, target, remaining - 1, tmpdir)
                    if dist == remaining - 1:
                        found = (move, candidate)
                        break
                if found is None:
                    return {
                        "status": "unavailable",
                        "message": f"could not certify step with {remaining} rotations remaining",
                        "steps": steps,
                    }
                move, current = found
                steps.append({"step": len(steps), "move": move, "tree": current})
                remaining -= 1

        return {"status": "ready", "steps": steps}

    def build_witnesses(self, rows: list[dict]) -> dict[str, dict]:
        witnesses: dict[str, dict] = {}
        for row in rows:
            direction = str(row.get("direction", ""))
            if row.get("status") != "ok" or direction not in {"a->b", "b->a"}:
                continue
            try:
                witnesses[direction] = self.build_witness(
                    str(row["tree_a"]),
                    str(row["tree_b"]),
                    int(row["distance"]),
                )
            except Exception as exc:
                witnesses[direction] = {"status": "unavailable", "message": str(exc)}
        return witnesses

    def handle_generate(self, raw_query: str) -> None:
        if not self.solver.exists():
            self.send_json(503, {
                "error": "Missing flipdist binary. Build first with: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j",
                "solver": str(self.solver),
            })
            return

        query = urllib.parse.parse_qs(raw_query)
        try:
            n = self.parse_int(query, "n", 12, 5, 60)
            seed = self.parse_int(query, "seed", 0, -2147483648, 2147483647)
            max_k_raw = query.get("max_k", [""])[0].strip()
            max_k = None
            if max_k_raw:
                max_k = self.parse_int(query, "max_k", 0, 1, 10000)
            else:
                # Li-Xia/STT diameter budget for n internal tree nodes.
                # Clamp to 1 so tiny generated examples still pass CLI validation.
                max_k = max(1, 2 * n - 6)

            case_kind = query.get("case", ["hard"])[0].strip().lower()
            case_map = {
                "easy": "simple",
                "simple": "simple",
                "hard": "random",
                "complex": "random",
                "random": "random",
            }
            if case_kind not in case_map:
                raise ValueError("case must be easy, hard, simple, complex, or random.")
            solver_case = case_map[case_kind]
        except ValueError as exc:
            self.send_json(400, {"error": str(exc)})
            return

        cmd = [
            str(self.solver),
            "--case",
            solver_case,
            "--n",
            str(n),
            "--count",
            "1",
            "--bfs-cap",
            "1",
        ]
        if solver_case == "random":
            cmd.extend(["--seed", str(seed)])
        if max_k is not None:
            cmd.extend(["--max-k", str(max_k)])

        try:
            proc = subprocess.run(
                cmd,
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                check=False,
                timeout=self.solver_timeout,
            )
        except subprocess.TimeoutExpired:
            self.send_json(504, {
                "error": f"flipdist timed out after {self.solver_timeout:g}s.",
                "command": cmd,
            })
            return

        if proc.returncode != 0:
            self.send_json(500, {
                "error": "flipdist returned a nonzero exit code.",
                "returncode": proc.returncode,
                "stderr": proc.stderr,
                "command": cmd,
            })
            return

        rows = [json.loads(line) for line in proc.stdout.splitlines() if line.strip()]
        witnesses = self.build_witnesses(rows)

        self.send_json(200, {
            "jsonl": proc.stdout,
            "stderr": proc.stderr,
            "command": cmd,
            "case": case_kind,
            "solver_case": solver_case,
            "max_k": max_k,
            "witnesses": witnesses,
            "witness_note": "Visualizer-only certificate reconstructed by repeated calls to the existing solver; normal flipdist CLI and solver logic are unchanged.",
        })


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Serve the FlipDist browser visualizer.")
    parser.add_argument("--host", default="127.0.0.1", help="Host interface to bind.")
    parser.add_argument("--port", type=int, default=8765, help="Preferred port.")
    parser.add_argument("--solver", default=None, help="Path to the flipdist binary.")
    parser.add_argument("--solver-timeout", type=float, default=30.0, help="Seconds to wait for generated cases.")
    parser.add_argument("--no-open", action="store_true", help="Do not open a browser window.")
    return parser.parse_args()


def port_is_free(host: str, port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            sock.bind((host, port))
        except OSError:
            return False
    return True


def choose_port(host: str, preferred: int) -> int:
    for port in range(preferred, preferred + 50):
        if port_is_free(host, port):
            return port
    raise RuntimeError(f"No free port found from {preferred} to {preferred + 49}.")


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parent
    repo_root = root.parents[1]
    solver = Path(args.solver).expanduser() if args.solver else repo_root / "build" / "flipdist"
    if not solver.is_absolute():
        solver = (repo_root / solver).resolve()
    port = choose_port(args.host, args.port)
    handler = functools.partial(
        VisualizerHandler,
        directory=str(root),
        repo_root=repo_root,
        solver=solver,
        solver_timeout=args.solver_timeout,
    )

    with ReusableThreadingServer((args.host, port), handler) as httpd:
        url = f"http://{args.host}:{port}/"
        print(f"Serving FlipDist visualizer at {url}", flush=True)
        print(f"Using solver: {solver}", flush=True)
        print("Witness paths are visualizer-only certificates built from existing solver calls; solver logic and public CLI are unchanged.", flush=True)
        print("Press Ctrl-C to stop.", flush=True)
        if not args.no_open:
            webbrowser.open(url)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nStopping visualizer server.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
