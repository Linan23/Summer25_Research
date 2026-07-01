/*
  Binary-tree utilities for the visualizer.
  The code parses FlipDist's P:...;I:... encoding, applies BST rotations, and
  serializes states for browser-side path reconstruction.
*/
export function parseTreeEncoding(encoding) {
  if (typeof encoding !== "string") {
    throw new Error("Tree encoding must be a string.");
  }

  const parts = new Map();
  for (const raw of encoding.split(";")) {
    const idx = raw.indexOf(":");
    if (idx < 0) continue;
    parts.set(raw.slice(0, idx).trim(), raw.slice(idx + 1).trim());
  }

  if (!parts.has("P") || !parts.has("I")) {
    throw new Error(`Invalid tree encoding: ${encoding}`);
  }

  const parseList = (text) => {
    if (text === "") return [];
    return text.split(",").map((item) => {
      const value = Number.parseInt(item.trim(), 10);
      if (!Number.isFinite(value)) {
        throw new Error(`Invalid node label "${item}" in ${encoding}`);
      }
      return value;
    });
  };

  const preorder = parseList(parts.get("P"));
  const inorder = parseList(parts.get("I"));
  if (preorder.length !== inorder.length) {
    throw new Error("Preorder and inorder lengths do not match.");
  }
  if (new Set(preorder).size !== preorder.length || new Set(inorder).size !== inorder.length) {
    throw new Error("Tree encodings must not contain duplicate node labels.");
  }
  const inSet = new Set(inorder);
  for (const value of preorder) {
    if (!inSet.has(value)) {
      throw new Error(`Node ${value} appears in preorder but not inorder.`);
    }
  }

  return { preorder, inorder };
}

export function treeFromEncoding(encoding) {
  const parsed = parseTreeEncoding(encoding);
  return buildTree(parsed.preorder, parsed.inorder);
}

export function buildTree(preorder, inorder) {
  const positions = new Map();
  inorder.forEach((value, index) => positions.set(value, index));
  const nodes = new Map();

  function build(preStart, preEnd, inStart, inEnd, parent) {
    if (preStart > preEnd) return null;
    const value = preorder[preStart];
    const inIndex = positions.get(value);
    if (inIndex === undefined || inIndex < inStart || inIndex > inEnd) {
      throw new Error("Preorder/inorder data does not describe one valid tree.");
    }

    const leftSize = inIndex - inStart;
    const node = { value, left: null, right: null, parent };
    nodes.set(value, node);
    node.left = build(preStart + 1, preStart + leftSize, inStart, inIndex - 1, value);
    node.right = build(preStart + leftSize + 1, preEnd, inIndex + 1, inEnd, value);
    return value;
  }

  const root = preorder.length ? build(0, preorder.length - 1, 0, inorder.length - 1, null) : null;
  return { root, nodes, inorder: [...inorder] };
}

export function cloneTree(state) {
  const nodes = new Map();
  for (const [value, node] of state.nodes) {
    nodes.set(value, {
      value,
      left: node.left,
      right: node.right,
      parent: node.parent
    });
  }
  return { root: state.root, nodes, inorder: [...state.inorder] };
}

export function nodeCount(state) {
  return state.nodes.size;
}

export function sameInorder(a, b) {
  if (a.inorder.length !== b.inorder.length) return false;
  return a.inorder.every((value, index) => value === b.inorder[index]);
}

export function preorderValues(state) {
  const out = [];
  function visit(value) {
    if (value === null || value === undefined) return;
    const node = state.nodes.get(value);
    if (!node) return;
    out.push(value);
    visit(node.left);
    visit(node.right);
  }
  visit(state.root);
  return out;
}

export function keyForTree(state) {
  return preorderValues(state).join(",");
}

export function treeToEncoding(state) {
  return `P:${preorderValues(state).join(",")};I:${state.inorder.join(",")}`;
}

export function rotateLeft(state, pivotValue) {
  const pivot = state.nodes.get(pivotValue);
  if (!pivot || pivot.right === null || pivot.right === undefined) return false;

  const childValue = pivot.right;
  const child = state.nodes.get(childValue);
  const betaValue = child.left;
  const parentValue = pivot.parent;

  child.parent = parentValue;
  if (parentValue === null || parentValue === undefined) {
    state.root = childValue;
  } else {
    const parent = state.nodes.get(parentValue);
    if (parent.left === pivotValue) parent.left = childValue;
    else if (parent.right === pivotValue) parent.right = childValue;
  }

  child.left = pivotValue;
  pivot.parent = childValue;
  pivot.right = betaValue;
  if (betaValue !== null && betaValue !== undefined) {
    state.nodes.get(betaValue).parent = pivotValue;
  }

  return true;
}

export function rotateRight(state, pivotValue) {
  const pivot = state.nodes.get(pivotValue);
  if (!pivot || pivot.left === null || pivot.left === undefined) return false;

  const childValue = pivot.left;
  const child = state.nodes.get(childValue);
  const betaValue = child.right;
  const parentValue = pivot.parent;

  child.parent = parentValue;
  if (parentValue === null || parentValue === undefined) {
    state.root = childValue;
  } else {
    const parent = state.nodes.get(parentValue);
    if (parent.left === pivotValue) parent.left = childValue;
    else if (parent.right === pivotValue) parent.right = childValue;
  }

  child.right = pivotValue;
  pivot.parent = childValue;
  pivot.left = betaValue;
  if (betaValue !== null && betaValue !== undefined) {
    state.nodes.get(betaValue).parent = pivotValue;
  }

  return true;
}

export function edgeSet(state) {
  const edges = new Set();
  for (const node of state.nodes.values()) {
    if (node.left !== null && node.left !== undefined) {
      edges.add(`${node.value}->${node.left}`);
    }
    if (node.right !== null && node.right !== undefined) {
      edges.add(`${node.value}->${node.right}`);
    }
  }
  return edges;
}

export function changedCurrentEdges(previousState, currentState) {
  if (!previousState || !currentState) return new Set();
  const before = edgeSet(previousState);
  const after = edgeSet(currentState);
  const changed = new Set();
  for (const edge of after) {
    if (!before.has(edge)) changed.add(edge);
  }
  return changed;
}

export function generateRotations(state) {
  const out = [];
  for (const value of state.inorder) {
    const node = state.nodes.get(value);
    if (!node) continue;

    if (node.right !== null && node.right !== undefined) {
      const next = cloneTree(state);
      if (rotateLeft(next, value)) {
        out.push({
          state: next,
          key: keyForTree(next),
          move: `rotateLeft(${value})`,
          pivot: value
        });
      }
    }

    if (node.left !== null && node.left !== undefined) {
      const next = cloneTree(state);
      if (rotateRight(next, value)) {
        out.push({
          state: next,
          key: keyForTree(next),
          move: `rotateRight(${value})`,
          pivot: value
        });
      }
    }
  }
  return out;
}
