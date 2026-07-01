/*
  JSON-lines parsing and case grouping for normal flipdist output.
  The parser intentionally consumes the existing CLI rows and does not require
  any path-specific output mode.
*/
export function parseJsonLines(text) {
  const rows = [];
  const errors = [];
  const lines = text.split(/\r?\n/);

  lines.forEach((line, index) => {
    const trimmed = line.trim();
    if (!trimmed) return;
    try {
      rows.push({ ...JSON.parse(trimmed), __line: index + 1 });
    } catch (error) {
      errors.push(`Line ${index + 1}: ${error.message}`);
    }
  });

  return { rows, errors };
}

function canonicalPairKey(a, b) {
  return [a, b].sort().join("||");
}

function caseKey(row) {
  return [
    row.case_type ?? "unknown",
    row.n ?? "n?",
    row.seed ?? "seed?",
    canonicalPairKey(row.tree_a, row.tree_b)
  ].join("|");
}

function inferGlobalTrees(rows) {
  const ab = rows.find((row) => row.direction === "a->b");
  if (ab) return { treeA: ab.tree_a, treeB: ab.tree_b };

  const ba = rows.find((row) => row.direction === "b->a");
  if (ba) return { treeA: ba.tree_b, treeB: ba.tree_a };

  const first = rows[0];
  return { treeA: first.tree_a, treeB: first.tree_b };
}

export function groupCases(rows) {
  const groups = new Map();

  rows.forEach((row) => {
    if (!row || typeof row.tree_a !== "string" || typeof row.tree_b !== "string") {
      return;
    }
    const key = caseKey(row);
    if (!groups.has(key)) {
      groups.set(key, []);
    }
    groups.get(key).push(row);
  });

  return [...groups.values()].map((caseRows, index) => {
    const first = caseRows[0];
    const directions = new Map();
    for (const row of caseRows) {
      const direction = row.direction || `row-${row.__line ?? directions.size + 1}`;
      if (!directions.has(direction)) {
        directions.set(direction, row);
      }
    }

    const globalTrees = inferGlobalTrees(caseRows);
    return {
      id: index,
      rows: caseRows,
      directions,
      treeA: globalTrees.treeA,
      treeB: globalTrees.treeB,
      label: formatCaseLabel(first, index)
    };
  });
}

export function formatCaseLabel(row, index) {
  const type = row.case_type ?? "case";
  const n = row.n ?? "?";
  const seed = row.seed ?? "?";
  return `${index + 1}. ${type} n=${n} seed=${seed}`;
}
