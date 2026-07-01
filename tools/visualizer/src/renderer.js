/*
  SVG tree renderer for the visualizer.
  It uses a standard top-down BST layout: inorder position defines the
  left-to-right coordinate and depth defines the vertical coordinate.
*/
function escapeText(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function layoutTree(state) {
  const depths = new Map();
  let maxDepth = 0;

  function walk(value, depth) {
    if (value === null || value === undefined) return null;
    const node = state.nodes.get(value);
    if (!node) return null;

    maxDepth = Math.max(maxDepth, depth);
    depths.set(value, depth);
    walk(node.left, depth + 1);
    walk(node.right, depth + 1);
  }

  walk(state.root, 0);
  if (!depths.size) {
    return { positions: new Map(), width: 420, height: 300, radius: 18, maxDepth: 0 };
  }

  const nodeCount = state.nodes.size;
  const radius = nodeCount <= 15 ? 18 : nodeCount <= 30 ? 15 : 12;
  const xGap = nodeCount <= 15 ? 74 : nodeCount <= 30 ? 58 : 46;
  const yGap = nodeCount <= 15 ? 76 : nodeCount <= 30 ? 66 : 56;
  const marginX = radius + 28;
  const marginY = radius + 26;
  const width = Math.max(360, marginX * 2 + Math.max(1, state.inorder.length - 1) * xGap);
  const height = Math.max(280, marginY * 2 + maxDepth * yGap);
  const positions = new Map();

  state.inorder.forEach((value, index) => {
    if (!depths.has(value)) return;
    positions.set(value, {
      x: marginX + index * xGap,
      y: marginY + depths.get(value) * yGap
    });
  });

  return { positions, width, height, radius, maxDepth };
}

export function renderEmpty(container, message) {
  container.innerHTML = `<div class="empty-tree">${escapeText(message)}</div>`;
}

export function renderTree(container, state, options = {}) {
  if (!state || !state.nodes || state.nodes.size === 0) {
    renderEmpty(container, "No tree loaded.");
    return;
  }

  const { positions, width, height, radius, maxDepth } = layoutTree(state);
  const highlightNode = options.highlightNode ?? null;
  const changedEdges = options.changedEdges ?? new Set();
  const edgeParts = [];
  const nodeParts = [];

  for (const node of state.nodes.values()) {
    const from = positions.get(node.value);
    for (const childValue of [node.left, node.right]) {
      if (childValue === null || childValue === undefined) continue;
      const to = positions.get(childValue);
      if (!from || !to) continue;
      const edgeKey = `${node.value}->${childValue}`;
      const className = changedEdges.has(edgeKey) ? "edge changed" : "edge";
      edgeParts.push(
        `<line class="${className}" x1="${from.x.toFixed(2)}" y1="${from.y.toFixed(2)}" x2="${to.x.toFixed(2)}" y2="${to.y.toFixed(2)}"></line>`
      );
    }
  }

  for (const value of state.inorder) {
    const pos = positions.get(value);
    if (!pos) continue;
    const className = value === highlightNode ? "node highlight" : "node";
    nodeParts.push(
      `<g class="${className}" transform="translate(${pos.x.toFixed(2)},${pos.y.toFixed(2)})">` +
      `<circle r="${radius}"></circle>` +
      `<text>${escapeText(value)}</text>` +
      `</g>`
    );
  }

  container.innerHTML =
    `<svg class="tree-svg" viewBox="0 0 ${width.toFixed(2)} ${height.toFixed(2)}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Binary tree with ${state.nodes.size} nodes and depth ${maxDepth}">` +
    edgeParts.join("") +
    nodeParts.join("") +
    `</svg>`;
}
