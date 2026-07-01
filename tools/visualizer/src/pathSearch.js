/*
  Bounded browser-side rotation-path search.
  This reconstructs one shortest path for visualization when n is small. The
  C++ solver remains the source of truth for the exact distance.
*/
import {
  generateRotations,
  keyForTree,
  nodeCount,
  sameInorder,
  treeFromEncoding
} from "./tree.js";

function yieldToBrowser() {
  return new Promise((resolve) => window.setTimeout(resolve, 0));
}

function buildPathFromKeys(keys, states) {
  return keys.map((key, index) => {
    const currentState = states.get(key);
    if (index === 0) {
      return { state: currentState, move: "start", pivot: null };
    }

    const previousState = states.get(keys[index - 1]);
    const move = generateRotations(previousState).find((next) => next.key === key);
    return {
      state: currentState,
      move: move?.move ?? "rotation",
      pivot: move?.pivot ?? null
    };
  });
}

function buildPathKeys(meetKey, frontParents, backParents) {
  const left = [];
  for (let key = meetKey; key !== null; key = frontParents.get(key).parentKey) {
    left.push(key);
  }
  left.reverse();

  const right = [];
  for (let key = backParents.get(meetKey)?.parentKey ?? null; key !== null; key = backParents.get(key).parentKey) {
    right.push(key);
  }

  return left.concat(right);
}

async function findShortestRotationPath(source, target, distance, stateCap, chunkSize) {
  const sourceKey = keyForTree(source);
  const targetKey = keyForTree(target);
  if (sourceKey === targetKey) {
    if (distance !== 0) {
      return {
        status: "mismatch",
        message: `The trees match but the solver distance is ${distance}.`
      };
    }
    return { status: "ready", steps: [{ state: source, move: "start", pivot: null }], visited: 1 };
  }

  let expanded = 0;

  let front = new Set([sourceKey]);
  let back = new Set([targetKey]);
  const frontParents = new Map([[sourceKey, { parentKey: null, depth: 0 }]]);
  const backParents = new Map([[targetKey, { parentKey: null, depth: 0 }]]);
  const frontStates = new Map([[sourceKey, source]]);
  const backStates = new Map([[targetKey, target]]);

  function combinedStateMap() {
    return new Map([...frontStates, ...backStates]);
  }

  while (front.size && back.size) {
    const frontDepth = frontParents.get(front.values().next().value).depth;
    const backDepth = backParents.get(back.values().next().value).depth;
    if (frontDepth + backDepth >= distance) break;

    const expandForward = front.size <= back.size;
    const currentFrontier = expandForward ? front : back;
    const ownParents = expandForward ? frontParents : backParents;
    const otherParents = expandForward ? backParents : frontParents;
    const ownStates = expandForward ? frontStates : backStates;
    const otherStates = expandForward ? backStates : frontStates;
    const nextFrontier = new Set();

    for (const key of currentFrontier) {
      const parent = ownParents.get(key);
      const state = ownStates.get(key);
      const nextDepth = parent.depth + 1;

      if (nextDepth + (expandForward ? backDepth : frontDepth) > distance) continue;

      for (const next of generateRotations(state)) {
        if (ownParents.has(next.key)) continue;

        ownParents.set(next.key, { parentKey: key, depth: nextDepth });
        ownStates.set(next.key, next.state);

        if (otherParents.has(next.key)) {
          const totalDepth = nextDepth + otherParents.get(next.key).depth;
          if (totalDepth !== distance) {
            continue;
          }
          if (!otherStates.has(next.key)) {
            otherStates.set(next.key, next.state);
          }
          const pathKeys = buildPathKeys(next.key, frontParents, backParents);
          const states = combinedStateMap();
          states.set(next.key, next.state);
          return {
            status: "ready",
            steps: buildPathFromKeys(pathKeys, states),
            visited: frontParents.size + backParents.size
          };
        }

        nextFrontier.add(next.key);

        if (frontParents.size + backParents.size >= stateCap) {
          return {
            status: "capped",
            message: `Path search stopped after ${(frontParents.size + backParents.size).toLocaleString()} states. The exact solver result is still shown.`
          };
        }
      }

      expanded += 1;
      if (expanded % chunkSize === 0) {
        await yieldToBrowser();
      }
    }

    if (expandForward) {
      front = nextFrontier;
    } else {
      back = nextFrontier;
    }
  }

  return {
    status: "not_found",
    message: `No path was reconstructed at depth ${distance} before the search frontier ended.`
  };
}

export async function findRotationPath(sourceEncoding, targetEncoding, distance, options = {}) {
  const parsedMaxNodes = Number(options.maxNodes ?? 12);
  const parsedStateCap = Number(options.stateCap ?? 150000);
  const parsedChunkSize = Number(options.chunkSize ?? 1200);
  const maxNodes = Number.isFinite(parsedMaxNodes) && parsedMaxNodes > 0 ? parsedMaxNodes : 12;
  const stateCap = Number.isFinite(parsedStateCap) && parsedStateCap > 0 ? parsedStateCap : 150000;
  const chunkSize = Number.isFinite(parsedChunkSize) && parsedChunkSize > 0 ? parsedChunkSize : 1200;

  if (!Number.isInteger(distance) || distance < 0) {
    return { status: "unavailable", message: "Path search requires a nonnegative solver distance." };
  }

  const source = treeFromEncoding(sourceEncoding);
  const target = treeFromEncoding(targetEncoding);
  const n = nodeCount(source);
  if (nodeCount(target) !== n || !sameInorder(source, target)) {
    return { status: "error", message: "Source and target trees do not have the same inorder node set." };
  }

  if (n > maxNodes) {
    return {
      status: "too_large",
      message: `Shortest-path animation is disabled for n=${n}. Current exact animation limit is n<=${maxNodes}. The solver distance shown above is still exact.`
    };
  }

  const shortest = await findShortestRotationPath(source, target, distance, stateCap, chunkSize);
  return shortest.status === "ready" ? { ...shortest, mode: "shortest" } : shortest;
}
