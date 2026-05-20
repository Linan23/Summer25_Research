import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

public final class TriangulationMetrics {
    private TriangulationMetrics() {}

    private static final String DEFAULT_PROGRAM = "triangulation_java";
    private static final String DEFAULT_CASE = "unknown";
    private static final String DEFAULT_DIRECTION = "a->b";

    public static final class Options {
        public String program = DEFAULT_PROGRAM;
        public String caseType = DEFAULT_CASE;
        public String direction = DEFAULT_DIRECTION;
        public long seed = -1L;
        public double timeLimitSeconds = 5.0;
        public long visitedCap = 5_000_000L;
        public long queueCap = 5_000_000L;

        public Options copy() {
            Options o = new Options();
            o.program = this.program;
            o.caseType = this.caseType;
            o.direction = this.direction;
            o.seed = this.seed;
            o.timeLimitSeconds = this.timeLimitSeconds;
            o.visitedCap = this.visitedCap;
            o.queueCap = this.queueCap;
            return o;
        }
    }

    public static final class Result {
        public String program;
        public String caseType;
        public long n;
        public long seed;
        public String direction;
        public long distance = -1L;
        public double timeMs = 0.0;
        public long expanded = 0L;
        public long enqueued = 0L;
        public long visited = 0L;
        public long maxQueue = 0L;
        public long duplicates = 0L;
        public String status = "ok";
        public String treeA;
        public String treeB;
        public String solver = "bfs";

        public String toJson() {
            StringBuilder sb = new StringBuilder();
            sb.append("{\"program\":\"").append(escapeJson(program)).append('"');
            sb.append(",\"case_type\":\"").append(escapeJson(caseType)).append('"');
            sb.append(",\"n\":").append(n);
            sb.append(",\"seed\":").append(seed);
            sb.append(",\"direction\":\"").append(escapeJson(direction)).append('"');
            sb.append(",\"distance\":").append(distance);
            sb.append(",\"time_ms\":").append(String.format(Locale.ROOT, "%.3f", timeMs));
            sb.append(",\"expanded\":").append(expanded);
            sb.append(",\"enqueued\":").append(enqueued);
            sb.append(",\"visited\":").append(visited);
            sb.append(",\"max_queue\":").append(maxQueue);
            sb.append(",\"duplicates\":").append(duplicates);
            sb.append(",\"status\":\"").append(escapeJson(status)).append('"');
            sb.append(",\"solver\":\"").append(escapeJson(solver)).append('"');
            sb.append(",\"tree_a\":\"").append(escapeJson(treeA)).append('"');
            sb.append(",\"tree_b\":\"").append(escapeJson(treeB)).append("\"}");
            return sb.toString();
        }
    }

    public static final class ResultPair {
        public final Result forward;
        public final Result reverse;
        public final boolean distanceAgrees;

        public ResultPair(Result forward, Result reverse, boolean distanceAgrees) {
            this.forward = forward;
            this.reverse = reverse;
            this.distanceAgrees = distanceAgrees;
        }
    }

    private enum StepOutcome { CONTINUE, FOUND, ABORT }

    public static Result run(Triangulation start, Triangulation target, Options options) {
        Options opts = (options == null) ? new Options() : options;

        Result res = new Result();
        res.program = safe(opts.program, DEFAULT_PROGRAM);
        res.caseType = safe(opts.caseType, DEFAULT_CASE);
        res.direction = safe(opts.direction, DEFAULT_DIRECTION);
        res.seed = opts.seed;

        Triangulation startCopy = start.clone2();
        Triangulation targetCopy = target.clone2();

        res.treeA = startCopy.canonicalTraversal();
        res.treeB = targetCopy.canonicalTraversal();
        res.n = startCopy.getInternalNodeIds().size();

        final String targetRep = targetCopy.getRep();
        final String startRep = startCopy.getRep();

        long startNs = System.nanoTime();
        boolean enforceTime = opts.timeLimitSeconds > 0.0 && Double.isFinite(opts.timeLimitSeconds);
        long deadline = enforceTime
            ? startNs + (long)Math.floor(opts.timeLimitSeconds * 1_000_000_000L)
            : Long.MAX_VALUE;

        if (startRep.equals(targetRep)) {
            res.distance = 0L;
            res.enqueued = 1L;
            res.visited = 1L;
            res.maxQueue = 1L;
            res.timeMs = (System.nanoTime() - startNs) / 1_000_000.0;
            return res;
        }

        ArrayDeque<Triangulation> queue = new ArrayDeque<>();
        Map<String, Integer> distMap = new HashMap<>();
        queue.add(startCopy);
        distMap.put(startRep, 0);
        res.enqueued = 1L;
        res.maxQueue = 1L;

        long visitedCap = (opts.visitedCap <= 0L) ? Long.MAX_VALUE : opts.visitedCap;
        long queueCap = (opts.queueCap <= 0L) ? Long.MAX_VALUE : opts.queueCap;

        boolean found = false;

        outer:
        while (!queue.isEmpty()) {
            if (enforceTime && System.nanoTime() >= deadline) {
                res.status = "timeout";
                break;
            }

            Triangulation cur = queue.poll();
            String curRep = cur.getRep();
            Integer depth = distMap.get(curRep);
            int curDist = (depth == null) ? 0 : depth;
            res.expanded++;

            List<Integer> ids = cur.getInternalNodeIds();
            for (int id : ids) {
                if (enforceTime && System.nanoTime() >= deadline) {
                    res.status = "timeout";
                    break outer;
                }

                if (cur.canRotateRightById(id)) {
                    Triangulation nxt = cur.clone2();
                    nxt.rotateRightById(id);
                    StepOutcome outcome = processNeighbor(nxt, curDist + 1, targetRep,
                                                          distMap, queue, res,
                                                          visitedCap, queueCap);
                    if (outcome == StepOutcome.FOUND) {
                        found = true;
                        break outer;
                    } else if (outcome == StepOutcome.ABORT) {
                        break outer;
                    }
                }
                if (cur.canRotateLeftById(id)) {
                    Triangulation nxt = cur.clone2();
                    nxt.rotateLeftById(id);
                    StepOutcome outcome = processNeighbor(nxt, curDist + 1, targetRep,
                                                          distMap, queue, res,
                                                          visitedCap, queueCap);
                    if (outcome == StepOutcome.FOUND) {
                        found = true;
                        break outer;
                    } else if (outcome == StepOutcome.ABORT) {
                        break outer;
                    }
                }
            }

            if ("timeout".equals(res.status) || "cap".equals(res.status)) {
                break;
            }
        }

        res.visited = distMap.size();
        res.timeMs = (System.nanoTime() - startNs) / 1_000_000.0;

        if (!found && "ok".equals(res.status)) {
            res.status = "error:not_found";
            res.distance = -1L;
        }

        return res;
    }

    public static ResultPair runBidirectional(Triangulation a,
                                              Triangulation b,
                                              Options options) {
        Options base = (options == null) ? new Options() : options;

        Options forwardOpts = base.copy();
        forwardOpts.direction = "a->b";
        Result forward = run(a, b, forwardOpts);

        Options reverseOpts = base.copy();
        reverseOpts.direction = "b->a";
        Result reverse = run(b, a, reverseOpts);

        boolean agree = "ok".equals(forward.status)
            && "ok".equals(reverse.status)
            && forward.distance >= 0
            && forward.distance == reverse.distance;

        return new ResultPair(forward, reverse, agree);
    }

    private static StepOutcome processNeighbor(Triangulation nxt,
                                               int nextDepth,
                                               String targetRep,
                                               Map<String, Integer> distMap,
                                               ArrayDeque<Triangulation> queue,
                                               Result res,
                                               long visitedCap,
                                               long queueCap) {
        String key = nxt.getRep();

        if (distMap.containsKey(key)) {
            res.duplicates++;
            return StepOutcome.CONTINUE;
        }

        if (key.equals(targetRep)) {
            distMap.put(key, nextDepth);
            res.visited = distMap.size();
            res.distance = nextDepth;
            return StepOutcome.FOUND;
        }

        if ((long)distMap.size() >= visitedCap) {
            res.status = "cap";
            return StepOutcome.ABORT;
        }
        if ((long)queue.size() >= queueCap) {
            res.status = "cap";
            return StepOutcome.ABORT;
        }

        distMap.put(key, nextDepth);
        res.visited = distMap.size();
        queue.add(nxt);
        res.enqueued++;
        long qSize = queue.size();
        if (qSize > res.maxQueue) res.maxQueue = qSize;

        return StepOutcome.CONTINUE;
    }

    private static String safe(String value, String fallback) {
        if (value == null || value.isEmpty()) return fallback;
        return value;
    }

    private static String escapeJson(String s) {
        if (s == null) return "";
        StringBuilder out = new StringBuilder(s.length() + 16);
        for (int i = 0; i < s.length(); ++i) {
            char c = s.charAt(i);
            switch (c) {
                case '"':  out.append("\\\""); break;
                case '\\': out.append("\\\\"); break;
                case '\b': out.append("\\b"); break;
                case '\f': out.append("\\f"); break;
                case '\n': out.append("\\n"); break;
                case '\r': out.append("\\r"); break;
                case '\t': out.append("\\t"); break;
                default:
                    if (c < 0x20) {
                        out.append("\\u");
                        String hex = Integer.toHexString(c);
                        for (int pad = hex.length(); pad < 4; pad++) out.append('0');
                        out.append(hex);
                    } else {
                        out.append(c);
                    }
            }
        }
        return out.toString();
    }
}
