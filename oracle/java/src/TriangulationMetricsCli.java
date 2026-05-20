import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

public class TriangulationMetricsCli {
    private static final double DEFAULT_TIME_LIMIT = 5.0;
    private static final long DEFAULT_VISITED_CAP = 5_000_000L;
    private static final long DEFAULT_QUEUE_CAP = 5_000_000L;

    public static void main(String[] args) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            try {
                Map<String, String> payload = parseJson(line);
                handleCase(payload);
            } catch (IllegalArgumentException ex) {
                System.err.println("Input parse error: " + ex.getMessage());
            }
        }
    }

    private static void handleCase(Map<String, String> payload) {
        String treeA = payload.get("tree_a");
        String treeB = payload.get("tree_b");
        if (treeA == null || treeB == null) {
            throw new IllegalArgumentException("Missing tree_a or tree_b");
        }

        String caseType = payload.getOrDefault("case_type", "unknown");
        long seed = payload.containsKey("seed") ? Long.parseLong(payload.get("seed")) : -1L;
        double timeLimit = payload.containsKey("time_limit")
            ? Double.parseDouble(payload.get("time_limit")) : DEFAULT_TIME_LIMIT;
        long visitedCap = payload.containsKey("visited_cap")
            ? Long.parseLong(payload.get("visited_cap")) : DEFAULT_VISITED_CAP;
        long queueCap = payload.containsKey("queue_cap")
            ? Long.parseLong(payload.get("queue_cap")) : DEFAULT_QUEUE_CAP;

        Triangulation start = Triangulation.fromCanonicalTraversal(treeA);
        Triangulation target = Triangulation.fromCanonicalTraversal(treeB);

        TriangulationMetrics.Options opts = new TriangulationMetrics.Options();
        opts.caseType = caseType;
        opts.seed = seed;
        opts.timeLimitSeconds = timeLimit;
        opts.visitedCap = visitedCap;
        opts.queueCap = queueCap;

        TriangulationMetrics.ResultPair pair = TriangulationMetrics.runBidirectional(start, target, opts);

        if (!pair.distanceAgrees &&
            pair.forward.distance >= 0 && pair.reverse.distance >= 0) {
            System.err.println(String.format(Locale.ROOT,
                "[WARN] distance mismatch seed=%d case=%s", seed, caseType));
        }

        System.out.println(pair.forward.toJson());
        System.out.println(pair.reverse.toJson());
    }

    private static Map<String, String> parseJson(String line) {
        String trimmed = line.trim();
        if (!trimmed.startsWith("{") || !trimmed.endsWith("}")) {
            throw new IllegalArgumentException("Expected JSON object: " + line);
        }
        String body = trimmed.substring(1, trimmed.length() - 1);
        List<String> tokens = splitTopLevel(body);

        Map<String, String> map = new HashMap<>();
        for (String token : tokens) {
            if (token.isEmpty()) continue;
            int colon = findColon(token);
            if (colon < 0) {
                throw new IllegalArgumentException("Missing ':' in token: " + token);
            }
            String rawKey = token.substring(0, colon).trim();
            String rawValue = token.substring(colon + 1).trim();

            if (!rawKey.startsWith("\"") || !rawKey.endsWith("\"")) {
                throw new IllegalArgumentException("Key must be quoted: " + rawKey);
            }
            String key = unescape(rawKey.substring(1, rawKey.length() - 1));
            String value;
            if (rawValue.startsWith("\"") && rawValue.endsWith("\"")) {
                value = unescape(rawValue.substring(1, rawValue.length() - 1));
            } else {
                value = rawValue;
            }
            map.put(key, value);
        }
        return map;
    }

    private static List<String> splitTopLevel(String body) {
        List<String> parts = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inString = false;
        boolean escaped = false;
        for (int i = 0; i < body.length(); ++i) {
            char c = body.charAt(i);
            if (escaped) {
                current.append(c);
                escaped = false;
                continue;
            }
            if (c == '\\') {
                current.append(c);
                escaped = true;
                continue;
            }
            if (c == '"') {
                current.append(c);
                inString = !inString;
                continue;
            }
            if (c == ',' && !inString) {
                parts.add(current.toString().trim());
                current.setLength(0);
                continue;
            }
            current.append(c);
        }
        if (current.length() > 0) {
            parts.add(current.toString().trim());
        }
        return parts;
    }

    private static int findColon(String token) {
        boolean inString = false;
        boolean escaped = false;
        for (int i = 0; i < token.length(); ++i) {
            char c = token.charAt(i);
            if (escaped) {
                escaped = false;
                continue;
            }
            if (c == '\\') {
                escaped = true;
                continue;
            }
            if (c == '"') {
                inString = !inString;
                continue;
            }
            if (c == ':' && !inString) {
                return i;
            }
        }
        return -1;
    }

    private static String unescape(String value) {
        StringBuilder out = new StringBuilder();
        boolean escaped = false;
        for (int i = 0; i < value.length(); ++i) {
            char c = value.charAt(i);
            if (!escaped) {
                if (c == '\\') {
                    escaped = true;
                } else {
                    out.append(c);
                }
                continue;
            }
            switch (c) {
                case '"': out.append('"'); break;
                case '\\': out.append('\\'); break;
                case 'n': out.append('\n'); break;
                case 'r': out.append('\r'); break;
                case 't': out.append('\t'); break;
                case 'b': out.append('\b'); break;
                case 'f': out.append('\f'); break;
                case 'u':
                    if (i + 4 >= value.length()) {
                        throw new IllegalArgumentException("Invalid unicode escape");
                    }
                    String hex = value.substring(i + 1, i + 5);
                    out.append((char) Integer.parseInt(hex, 16));
                    i += 4;
                    break;
                default:
                    out.append(c);
                    break;
            }
            escaped = false;
        }
        if (escaped) {
            throw new IllegalArgumentException("Trailing backslash in string");
        }
        return out.toString();
    }
}
