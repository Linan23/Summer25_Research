import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

public final class TriangulationFileDistanceCli {
    private static final double DEFAULT_TIME_LIMIT_SECONDS = 30.0;
    private static final long DEFAULT_VISITED_CAP = 2_000_000L;
    private static final long DEFAULT_QUEUE_CAP = 2_000_000L;

    private TriangulationFileDistanceCli() {}

    private static final class Point {
        final double x;
        final double y;

        Point(double x, double y) {
            this.x = x;
            this.y = y;
        }
    }

    private static final class EdgeKey {
        final int a;
        final int b;

        EdgeKey(int u, int v) {
            if (u <= v) {
                this.a = u;
                this.b = v;
            } else {
                this.a = v;
                this.b = u;
            }
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) return true;
            if (!(other instanceof EdgeKey)) return false;
            EdgeKey o = (EdgeKey) other;
            return a == o.a && b == o.b;
        }

        @Override
        public int hashCode() {
            return 31 * a + b;
        }
    }

    private static final class State {
        final int[][] triangles;
        final String key;

        State(int[][] triangles) {
            this.triangles = canonicalizeTriangles(triangles);
            this.key = encodeTriangles(this.triangles);
        }
    }

    private static final class SearchResult {
        long distance = -1L;
        String status = "ok";
        double timeMs = 0.0;
        long expanded = 0L;
        long visited = 0L;
        long enqueued = 0L;
        long maxQueue = 0L;

        String toJson(String pointsPath, String tri1Path, String tri2Path) {
            return String.format(Locale.ROOT,
                "{\"program\":\"triangulation_java_file_bfs\",\"solver\":\"bidir_bfs\"," +
                "\"points_path\":\"%s\",\"triangulation_1_path\":\"%s\",\"triangulation_2_path\":\"%s\"," +
                "\"distance\":%d,\"status\":\"%s\",\"time_ms\":%.3f,\"expanded\":%d," +
                "\"visited\":%d,\"enqueued\":%d,\"max_queue\":%d}",
                escapeJson(pointsPath),
                escapeJson(tri1Path),
                escapeJson(tri2Path),
                distance,
                escapeJson(status),
                timeMs,
                expanded,
                visited,
                enqueued,
                maxQueue
            );
        }
    }

    private static final class Frontier {
        final ArrayDeque<State> queue = new ArrayDeque<>();
        final Map<String, Integer> depth = new HashMap<>();

        Frontier(State start) {
            queue.add(start);
            depth.put(start.key, 0);
        }
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 3) {
            System.err.println("Usage: TriangulationFileDistanceCli points_path triangulation_1_path triangulation_2_path " +
                               "[--time-limit-sec N] [--visited-cap N] [--queue-cap N]");
            System.exit(2);
        }

        String pointsPath = args[0];
        String tri1Path = args[1];
        String tri2Path = args[2];
        double timeLimitSeconds = DEFAULT_TIME_LIMIT_SECONDS;
        long visitedCap = DEFAULT_VISITED_CAP;
        long queueCap = DEFAULT_QUEUE_CAP;

        for (int i = 3; i < args.length; ++i) {
            String flag = args[i];
            if ("--time-limit-sec".equals(flag) && i + 1 < args.length) {
                timeLimitSeconds = Double.parseDouble(args[++i]);
            } else if ("--visited-cap".equals(flag) && i + 1 < args.length) {
                visitedCap = Long.parseLong(args[++i]);
            } else if ("--queue-cap".equals(flag) && i + 1 < args.length) {
                queueCap = Long.parseLong(args[++i]);
            } else {
                throw new IllegalArgumentException("Unknown or incomplete flag: " + flag);
            }
        }

        Point[] points = loadPoints(Path.of(pointsPath));
        State start = new State(loadTriangles(Path.of(tri1Path)));
        State target = new State(loadTriangles(Path.of(tri2Path)));

        SearchResult result = bidirectionalDistance(points, start, target, timeLimitSeconds, visitedCap, queueCap);
        System.out.println(result.toJson(pointsPath, tri1Path, tri2Path));
    }

    private static SearchResult bidirectionalDistance(Point[] points,
                                                      State start,
                                                      State target,
                                                      double timeLimitSeconds,
                                                      long visitedCap,
                                                      long queueCap) {
        SearchResult result = new SearchResult();
        long startNs = System.nanoTime();
        boolean enforceTime = timeLimitSeconds > 0.0 && Double.isFinite(timeLimitSeconds);
        long deadlineNs = enforceTime
            ? startNs + (long) Math.floor(timeLimitSeconds * 1_000_000_000L)
            : Long.MAX_VALUE;

        if (start.key.equals(target.key)) {
            result.distance = 0L;
            result.visited = 2L;
            result.enqueued = 2L;
            result.maxQueue = 1L;
            result.timeMs = elapsedMs(startNs);
            return result;
        }

        Frontier forward = new Frontier(start);
        Frontier reverse = new Frontier(target);
        long maxVisited = (visitedCap <= 0L) ? Long.MAX_VALUE : visitedCap;
        long maxQueue = (queueCap <= 0L) ? Long.MAX_VALUE : queueCap;
        result.enqueued = 2L;
        result.visited = 2L;
        result.maxQueue = 1L;

        while (!forward.queue.isEmpty() && !reverse.queue.isEmpty()) {
            if (enforceTime && System.nanoTime() >= deadlineNs) {
                result.status = "timeout";
                break;
            }
            boolean expandForward = forward.queue.size() <= reverse.queue.size();
            int found = expandLayer(points,
                                    expandForward ? forward : reverse,
                                    expandForward ? reverse : forward,
                                    result,
                                    maxVisited,
                                    maxQueue,
                                    deadlineNs,
                                    enforceTime);
            if (found >= 0) {
                result.distance = found;
                break;
            }
            if (!"ok".equals(result.status)) {
                break;
            }
        }

        if (result.distance < 0 && "ok".equals(result.status)) {
            result.status = "error:not_found";
        }
        result.timeMs = elapsedMs(startNs);
        return result;
    }

    private static int expandLayer(Point[] points,
                                   Frontier active,
                                   Frontier other,
                                   SearchResult result,
                                   long visitedCap,
                                   long queueCap,
                                   long deadlineNs,
                                   boolean enforceTime) {
        int layerSize = active.queue.size();
        for (int i = 0; i < layerSize; ++i) {
            if (enforceTime && System.nanoTime() >= deadlineNs) {
                result.status = "timeout";
                return -1;
            }

            State cur = active.queue.poll();
            if (cur == null) {
                continue;
            }
            int curDepth = active.depth.getOrDefault(cur.key, 0);
            result.expanded++;

            for (State nxt : neighbors(points, cur)) {
                if (active.depth.containsKey(nxt.key)) {
                    continue;
                }
                Integer otherDepth = other.depth.get(nxt.key);
                if (otherDepth != null) {
                    return curDepth + 1 + otherDepth;
                }

                long combinedVisited = (long) active.depth.size() + (long) other.depth.size();
                if (combinedVisited >= visitedCap) {
                    result.status = "cap";
                    return -1;
                }
                long combinedQueue = (long) active.queue.size() + (long) other.queue.size();
                if (combinedQueue >= queueCap) {
                    result.status = "cap";
                    return -1;
                }

                active.depth.put(nxt.key, curDepth + 1);
                active.queue.add(nxt);
                result.visited = (long) active.depth.size() + (long) other.depth.size();
                result.enqueued++;
                long qSize = (long) active.queue.size() + (long) other.queue.size();
                if (qSize > result.maxQueue) {
                    result.maxQueue = qSize;
                }
            }
        }
        return -1;
    }

    private static List<State> neighbors(Point[] points, State state) {
        Map<EdgeKey, List<Integer>> edgeToTriangles = new HashMap<>();
        Set<EdgeKey> presentEdges = new HashSet<>();

        for (int i = 0; i < state.triangles.length; ++i) {
            int[] t = state.triangles[i];
            addEdge(edgeToTriangles, presentEdges, t[0], t[1], i);
            addEdge(edgeToTriangles, presentEdges, t[1], t[2], i);
            addEdge(edgeToTriangles, presentEdges, t[0], t[2], i);
        }

        List<State> out = new ArrayList<>();
        Set<String> seen = new HashSet<>();
        for (Map.Entry<EdgeKey, List<Integer>> entry : edgeToTriangles.entrySet()) {
            List<Integer> owners = entry.getValue();
            if (owners.size() != 2) {
                continue;
            }

            int idx1 = owners.get(0);
            int idx2 = owners.get(1);
            int[] t1 = state.triangles[idx1];
            int[] t2 = state.triangles[idx2];
            int u = entry.getKey().a;
            int v = entry.getKey().b;
            int c = thirdVertex(t1, u, v);
            int d = thirdVertex(t2, u, v);
            if (c < 0 || d < 0 || c == d) {
                continue;
            }
            EdgeKey newDiag = new EdgeKey(c, d);
            if (presentEdges.contains(newDiag)) {
                continue;
            }
            if (!segmentsProperlyIntersect(points[u], points[v], points[c], points[d])) {
                continue;
            }

            int[][] nextTriangles = new int[state.triangles.length][3];
            int write = 0;
            for (int i = 0; i < state.triangles.length; ++i) {
                if (i == idx1 || i == idx2) {
                    continue;
                }
                nextTriangles[write++] = state.triangles[i].clone();
            }
            nextTriangles[write++] = sortTriple(c, d, u);
            nextTriangles[write] = sortTriple(c, d, v);

            State nextState = new State(nextTriangles);
            if (seen.add(nextState.key)) {
                out.add(nextState);
            }
        }
        return out;
    }

    private static void addEdge(Map<EdgeKey, List<Integer>> edgeToTriangles,
                                Set<EdgeKey> presentEdges,
                                int u,
                                int v,
                                int triangleIndex) {
        EdgeKey key = new EdgeKey(u, v);
        presentEdges.add(key);
        edgeToTriangles.computeIfAbsent(key, ignored -> new ArrayList<>(2)).add(triangleIndex);
    }

    private static int thirdVertex(int[] triangle, int u, int v) {
        for (int value : triangle) {
            if (value != u && value != v) {
                return value;
            }
        }
        return -1;
    }

    private static boolean segmentsProperlyIntersect(Point a,
                                                     Point b,
                                                     Point c,
                                                     Point d) {
        double abC = orientation(a, b, c);
        double abD = orientation(a, b, d);
        double cdA = orientation(c, d, a);
        double cdB = orientation(c, d, b);
        return sign(abC) * sign(abD) < 0 && sign(cdA) * sign(cdB) < 0;
    }

    private static double orientation(Point a, Point b, Point c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    private static int sign(double value) {
        if (value > 0.0) return 1;
        if (value < 0.0) return -1;
        return 0;
    }

    private static Point[] loadPoints(Path path) throws IOException {
        List<Point> points = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\\s+");
                if (parts.length < 2) {
                    throw new IOException("Invalid points line: " + line);
                }
                points.add(new Point(Double.parseDouble(parts[0]), Double.parseDouble(parts[1])));
            }
        }
        return points.toArray(new Point[0]);
    }

    private static int[][] loadTriangles(Path path) throws IOException {
        List<int[]> triangles = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\\s+");
                if (parts.length < 3) {
                    throw new IOException("Invalid triangle line: " + line);
                }
                int a = Integer.parseInt(parts[0]);
                int b = Integer.parseInt(parts[1]);
                int c = Integer.parseInt(parts[2]);
                triangles.add(sortTriple(a, b, c));
            }
        }
        return triangles.toArray(new int[0][3]);
    }

    private static int[] sortTriple(int a, int b, int c) {
        int[] triple = new int[] { a, b, c };
        Arrays.sort(triple);
        return triple;
    }

    private static int[][] canonicalizeTriangles(int[][] triangles) {
        int[][] copy = new int[triangles.length][3];
        for (int i = 0; i < triangles.length; ++i) {
            copy[i] = triangles[i].clone();
            Arrays.sort(copy[i]);
        }
        Arrays.sort(copy, (left, right) -> {
            if (left[0] != right[0]) return Integer.compare(left[0], right[0]);
            if (left[1] != right[1]) return Integer.compare(left[1], right[1]);
            return Integer.compare(left[2], right[2]);
        });
        return copy;
    }

    private static String encodeTriangles(int[][] triangles) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < triangles.length; ++i) {
            if (i > 0) {
                sb.append(';');
            }
            sb.append(triangles[i][0]).append(',')
              .append(triangles[i][1]).append(',')
              .append(triangles[i][2]);
        }
        return sb.toString();
    }

    private static double elapsedMs(long startNs) {
        return (System.nanoTime() - startNs) / 1_000_000.0;
    }

    private static String escapeJson(String value) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < value.length(); ++i) {
            char c = value.charAt(i);
            switch (c) {
                case '\\':
                    sb.append("\\\\");
                    break;
                case '"':
                    sb.append("\\\"");
                    break;
                case '\n':
                    sb.append("\\n");
                    break;
                case '\r':
                    sb.append("\\r");
                    break;
                case '\t':
                    sb.append("\\t");
                    break;
                default:
                    sb.append(c);
                    break;
            }
        }
        return sb.toString();
    }
}
