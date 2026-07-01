/*
  Main browser application for the FlipDist visualizer.
  The app asks the local server to generate flipdist JSON-lines rows, shows
  both directions, and uses server witnesses or browser BFS for path animation.
*/
import { findRotationPath } from "./pathSearch.js";
import { groupCases, parseJsonLines } from "./parser.js";
import { renderEmpty, renderTree } from "./renderer.js";
import { TERMS } from "./terminology.js";
import {
  changedCurrentEdges,
  treeFromEncoding
} from "./tree.js";

const els = {
  generateN: document.getElementById("generateN"),
  generateCaseType: document.getElementById("generateCaseType"),
  generateSeed: document.getElementById("generateSeed"),
  generateMaxK: document.getElementById("generateMaxK"),
  generateCase: document.getElementById("generateCase"),
  caseSelect: document.getElementById("caseSelect"),
  dirAB: document.getElementById("dirAB"),
  dirBA: document.getElementById("dirBA"),
  stateCap: document.getElementById("stateCap"),
  caseType: document.getElementById("caseType"),
  caseN: document.getElementById("caseN"),
  caseSeed: document.getElementById("caseSeed"),
  caseDirection: document.getElementById("caseDirection"),
  caseDistance: document.getElementById("caseDistance"),
  caseStatus: document.getElementById("caseStatus"),
  caseRuntime: document.getElementById("caseRuntime"),
  caseStep: document.getElementById("caseStep"),
  messageLine: document.getElementById("messageLine"),
  playbackBand: document.getElementById("playbackBand"),
  resetStep: document.getElementById("resetStep"),
  prevStep: document.getElementById("prevStep"),
  playPause: document.getElementById("playPause"),
  nextStep: document.getElementById("nextStep"),
  stepSlider: document.getElementById("stepSlider"),
  moveLabel: document.getElementById("moveLabel"),
  treeA: document.getElementById("treeA"),
  treeB: document.getElementById("treeB"),
  treeCurrent: document.getElementById("treeCurrent"),
  treeGrid: document.getElementById("treeGrid"),
  terminologyList: document.getElementById("terminologyList")
};

const state = {
  cases: [],
  selectedCaseIndex: 0,
  selectedDirection: "a->b",
  pathCache: new Map(),
  witnesses: {},
  currentStep: 0,
  playTimer: null,
  searchToken: 0,
  generating: false
};

function setMessage(text, level = "info") {
  els.messageLine.textContent = text;
  els.messageLine.classList.toggle("warning", level === "warning");
  els.messageLine.classList.toggle("error", level === "error");
}

function selectedCase() {
  return state.cases[state.selectedCaseIndex] ?? null;
}

function selectedRow() {
  const current = selectedCase();
  if (!current) return null;
  return current.directions.get(state.selectedDirection) ?? null;
}

function cacheKey(row) {
  return [
    row.case_type ?? "",
    row.n ?? "",
    row.seed ?? "",
    row.direction ?? "",
    row.tree_a,
    row.tree_b,
    row.distance ?? ""
  ].join("|");
}

function stopPlayback() {
  if (state.playTimer !== null) {
    window.clearInterval(state.playTimer);
    state.playTimer = null;
  }
  els.playPause.textContent = "Play";
}

function renderTerms() {
  els.terminologyList.innerHTML = TERMS.map((entry) => (
    `<div><dt>${entry.term}</dt><dd>${entry.text}</dd></div>`
  )).join("");
}

function setMetricValues(row) {
  els.caseType.textContent = row?.case_type ?? "-";
  els.caseN.textContent = row?.n ?? "-";
  els.caseSeed.textContent = row?.seed ?? "-";
  els.caseDirection.textContent = row?.direction ?? state.selectedDirection ?? "-";
  els.caseDistance.textContent = row?.distance ?? "-";
  els.caseStatus.textContent = row?.status ?? "-";
  els.caseRuntime.textContent = row?.time_ms !== undefined ? `${row.time_ms} ms` : "-";
}

function renderReferenceTrees(currentCase) {
  if (!currentCase) {
    renderEmpty(els.treeA, "No case loaded.");
    renderEmpty(els.treeB, "No case loaded.");
    return;
  }

  try {
    renderTree(els.treeA, treeFromEncoding(currentCase.treeA));
    renderTree(els.treeB, treeFromEncoding(currentCase.treeB));
  } catch (error) {
    renderEmpty(els.treeA, error.message);
    renderEmpty(els.treeB, error.message);
  }
}

function renderCurrentTree(pathInfo) {
  const row = selectedRow();
  if (!row) {
    renderEmpty(els.treeCurrent, "No direction row loaded.");
    return;
  }

  if (pathInfo?.status === "ready") {
    const steps = pathInfo.steps;
    state.currentStep = Math.min(state.currentStep, steps.length - 1);
    const current = steps[state.currentStep];
    const previous = state.currentStep > 0 ? steps[state.currentStep - 1].state : null;
    renderTree(els.treeCurrent, current.state, {
      highlightNode: current.pivot,
      changedEdges: changedCurrentEdges(previous, current.state)
    });
    els.stepSlider.max = String(Math.max(0, steps.length - 1));
    els.stepSlider.value = String(state.currentStep);
    els.caseStep.textContent = `${state.currentStep}/${steps.length - 1}`;
    els.moveLabel.textContent = `Move: ${current.move}`;
    return;
  }

  try {
    renderTree(els.treeCurrent, treeFromEncoding(row.tree_a));
  } catch (error) {
    renderEmpty(els.treeCurrent, error.message);
  }
  els.stepSlider.max = "0";
  els.stepSlider.value = "0";
  els.caseStep.textContent = "-";
  els.moveLabel.textContent = "Move: -";
}

function setPlaybackEnabled(enabled) {
  els.playbackBand.hidden = !enabled;
  els.resetStep.disabled = !enabled;
  els.prevStep.disabled = !enabled;
  els.playPause.disabled = !enabled;
  els.nextStep.disabled = !enabled;
  els.stepSlider.disabled = !enabled;
}

function setAnimationVisible(visible) {
  els.treeGrid.classList.toggle("no-current", !visible);
  if (!visible) {
    renderEmpty(els.treeCurrent, "No animation available.");
  }
}

function updateDirectionTabs(currentCase) {
  const hasAB = currentCase?.directions.has("a->b") ?? false;
  const hasBA = currentCase?.directions.has("b->a") ?? false;
  els.dirAB.disabled = !hasAB;
  els.dirBA.disabled = !hasBA;
  els.dirAB.classList.toggle("is-active", state.selectedDirection === "a->b");
  els.dirBA.classList.toggle("is-active", state.selectedDirection === "b->a");
}

async function ensurePath(row) {
  if (!row) return { status: "missing", message: "No selected direction row." };
  const key = cacheKey(row);
  if (state.pathCache.has(key)) return state.pathCache.get(key);

  if (row.status !== "ok") {
    const unavailable = {
      status: "unavailable",
      message: `Path animation requires status=ok. This row has status=${row.status}.`
    };
    state.pathCache.set(key, unavailable);
    return unavailable;
  }

  const distance = Number(row.distance);
  if (!Number.isInteger(distance) || distance < 0) {
    const unavailable = {
      status: "unavailable",
      message: "Path animation requires a nonnegative integer distance."
    };
    state.pathCache.set(key, unavailable);
    return unavailable;
  }

  const witness = state.witnesses?.[row.direction];
  if (witness?.status === "ready" && Array.isArray(witness.steps)) {
    if (witness.steps.length === distance + 1) {
      const steps = witness.steps.map((step) => {
        const match = String(step.move).match(/\(([-0-9]+)\)/);
        return {
          state: treeFromEncoding(step.tree),
          move: step.move,
          pivot: match ? Number.parseInt(match[1], 10) : null
        };
      });
      const ready = {
        status: "ready",
        mode: "solver-witness",
        steps,
        visited: witness.steps.length,
        note: witness.note
      };
      state.pathCache.set(key, ready);
      return ready;
    }
  }

  const token = ++state.searchToken;
  const stateCap = Number.parseInt(els.stateCap.value, 10);
  setMessage("Reconstructing a visual rotation path in the browser...");
  const pathInfo = await findRotationPath(row.tree_a, row.tree_b, distance, {
    maxNodes: 14,
    stateCap: Number.isFinite(stateCap) && stateCap > 0 ? stateCap : 150000
  });

  if (token === state.searchToken) {
    state.pathCache.set(key, pathInfo);
  }
  return pathInfo;
}

function applyPathMessage(pathInfo) {
  const row = selectedRow();
  if (!pathInfo) {
    setMessage("No path information available.", "warning");
    return;
  }
  if (pathInfo.status === "ready") {
    if (pathInfo.mode === "solver-witness") {
      setMessage(`Shortest solver path ready from server witness: ${pathInfo.steps.length - 1} rotations.`);
      return;
    }
    setMessage(`Shortest solver path ready. Browser BFS visited ${pathInfo.visited.toLocaleString()} states.`);
    return;
  }
  if (row?.status === "ok") {
    setMessage("Solved exactly. Shortest-path animation is unavailable for this size.");
    return;
  }
  const level = pathInfo.status === "error" || pathInfo.status === "mismatch" ? "error" : "warning";
  setMessage(pathInfo.message ?? "Path animation unavailable for this case.", level);
}

async function renderAll() {
  stopPlayback();
  const currentCase = selectedCase();
  const row = selectedRow();

  renderReferenceTrees(currentCase);
  updateDirectionTabs(currentCase);
  setMetricValues(row);

  if (!currentCase || !row) {
    setPlaybackEnabled(false);
    setAnimationVisible(false);
    renderCurrentTree(null);
    setMessage("Choose n and a case type, then generate a solver case.", "warning");
    return;
  }

  state.currentStep = 0;
  setPlaybackEnabled(false);
  setAnimationVisible(false);
  renderCurrentTree(null);

  const pathInfo = await ensurePath(row);
  if (row !== selectedRow()) {
    return;
  }
  const hasAnimation = pathInfo.status === "ready";
  setAnimationVisible(hasAnimation);
  renderCurrentTree(hasAnimation ? pathInfo : null);
  setPlaybackEnabled(hasAnimation && pathInfo.steps.length > 1);
  applyPathMessage(pathInfo);
}

function populateCases(cases) {
  state.cases = cases;
  state.selectedCaseIndex = 0;
  state.selectedDirection = "a->b";
  state.pathCache.clear();
  els.caseSelect.innerHTML = "";

  cases.forEach((item, index) => {
    const option = document.createElement("option");
    option.value = String(index);
    option.textContent = item.label;
    els.caseSelect.appendChild(option);
  });

  const first = selectedCase();
  if (first && !first.directions.has("a->b")) {
    state.selectedDirection = first.directions.has("b->a") ? "b->a" : [...first.directions.keys()][0];
  }

  renderAll();
}

function useGeneratedOutput(label, text, witnesses = {}) {
  stopPlayback();
  state.witnesses = witnesses || {};
  const parsed = parseJsonLines(text);
  const cases = groupCases(parsed.rows);
  if (!cases.length) {
    setMessage(`No usable flipdist rows found for ${label}.`, "error");
    populateCases([]);
    return;
  }
  populateCases(cases);
  if (parsed.errors.length) {
    setMessage(`Generated ${cases.length} case(s) for ${label}; ignored ${parsed.errors.length} malformed line(s).`, "warning");
  } else {
    setMessage(`Generated ${cases.length} case(s) for ${label}.`);
  }
}

function updateGenerateCaseUi() {
  const isEasy = els.generateCaseType.value === "easy";
  els.generateSeed.disabled = isEasy;
  els.generateSeed.title = isEasy ? "Easy cases use the solver's simple-case generator, which does not use a seed." : "";
}

async function generateSolverCase() {
  if (state.generating) return;

  const n = Number.parseInt(els.generateN.value, 10);
  const seed = Number.parseInt(els.generateSeed.value, 10);
  const maxKText = els.generateMaxK.value.trim();

  if (!Number.isInteger(n) || n < 5 || n > 60) {
    setMessage("n must be an integer between 5 and 60.", "error");
    return;
  }
  if (!Number.isInteger(seed)) {
    setMessage("Seed must be an integer.", "error");
    return;
  }

  const params = new URLSearchParams({
    n: String(n),
    case: els.generateCaseType.value,
    seed: String(seed)
  });
  if (maxKText !== "") {
    const maxK = Number.parseInt(maxKText, 10);
    if (!Number.isInteger(maxK) || maxK < 1 || maxK > 10000) {
      setMessage("max_k must be blank or an integer between 1 and 10000.", "error");
      return;
    }
    params.set("max_k", String(maxK));
  }

  state.generating = true;
  els.generateCase.disabled = true;
  els.generateCase.textContent = "Generating";
  setMessage("Running flipdist locally through the visualizer server...");

  try {
    const response = await fetch(`api/generate?${params.toString()}`);
    let payload = {};
    try {
      payload = await response.json();
    } catch {
      payload = {};
    }

    if (!response.ok) {
      throw new Error(payload.error ?? `Generation failed with HTTP ${response.status}.`);
    }
    if (!payload.jsonl || !payload.jsonl.trim()) {
      throw new Error("flipdist returned no JSON-lines output.");
    }

    const seedLabel = els.generateCaseType.value === "easy" ? "" : ` seed=${seed}`;
    useGeneratedOutput(`${els.generateCaseType.value} n=${n}${seedLabel}`, payload.jsonl, payload.witnesses);
  } catch (error) {
    setMessage(`${error.message} Start with python3 tools/visualizer/serve.py if this page was opened directly.`, "error");
  } finally {
    state.generating = false;
    els.generateCase.disabled = false;
    els.generateCase.textContent = "Generate";
  }
}

function stepTo(index) {
  const row = selectedRow();
  if (!row) return;
  const pathInfo = state.pathCache.get(cacheKey(row));
  if (pathInfo?.status !== "ready") return;
  state.currentStep = Math.max(0, Math.min(index, pathInfo.steps.length - 1));
  renderCurrentTree(pathInfo);
}

function playOrPause() {
  const row = selectedRow();
  if (!row) return;
  const pathInfo = state.pathCache.get(cacheKey(row));
  if (pathInfo?.status !== "ready") return;

  if (state.playTimer !== null) {
    stopPlayback();
    return;
  }

  els.playPause.textContent = "Pause";
  state.playTimer = window.setInterval(() => {
    if (state.currentStep >= pathInfo.steps.length - 1) {
      stopPlayback();
      return;
    }
    stepTo(state.currentStep + 1);
  }, 650);
}

function bindEvents() {
  els.generateCase.addEventListener("click", generateSolverCase);
  els.generateCaseType.addEventListener("change", updateGenerateCaseUi);

  for (const input of [els.generateN, els.generateSeed, els.generateMaxK]) {
    input.addEventListener("keydown", (event) => {
      if (event.key === "Enter") {
        generateSolverCase();
      }
    });
  }

  els.caseSelect.addEventListener("change", () => {
    state.selectedCaseIndex = Number(els.caseSelect.value);
    const current = selectedCase();
    if (current && !current.directions.has(state.selectedDirection)) {
      state.selectedDirection = current.directions.has("a->b") ? "a->b" : [...current.directions.keys()][0];
    }
    renderAll();
  });

  for (const button of [els.dirAB, els.dirBA]) {
    button.addEventListener("click", () => {
      if (button.disabled) return;
      state.selectedDirection = button.dataset.direction;
      renderAll();
    });
  }

  els.stateCap.addEventListener("change", () => {
    state.pathCache.clear();
    renderAll();
  });

  els.resetStep.addEventListener("click", () => stepTo(0));
  els.prevStep.addEventListener("click", () => stepTo(state.currentStep - 1));
  els.nextStep.addEventListener("click", () => stepTo(state.currentStep + 1));
  els.playPause.addEventListener("click", playOrPause);
  els.stepSlider.addEventListener("input", () => stepTo(Number(els.stepSlider.value)));
}

function boot() {
  renderTerms();
  bindEvents();
  updateGenerateCaseUi();
  populateCases([]);
}

boot();
