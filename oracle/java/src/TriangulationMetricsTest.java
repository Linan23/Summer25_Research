import java.util.List;
import java.util.Random;

public class TriangulationMetricsTest {
    private static void assertTrue(boolean condition, String message) {
        if (!condition) {
            throw new AssertionError(message);
        }
    }

    private static void testIdentity() {
        Triangulation base = new Triangulation(4, new Random(42));
        Triangulation copy = base.clone2();

        TriangulationMetrics.Options opts = new TriangulationMetrics.Options();
        opts.timeLimitSeconds = 1.0;
        opts.visitedCap = 100_000;
        opts.queueCap = 100_000;

        TriangulationMetrics.ResultPair pair = TriangulationMetrics.runBidirectional(base, copy, opts);
        assertTrue(pair.distanceAgrees, "Identity pair should agree");
        assertTrue(pair.forward.distance == 0, "Identity distance must be zero");
        assertTrue("ok".equals(pair.forward.status), "Identity status should be ok");
        assertTrue("bfs".equals(pair.forward.solver), "Identity solver should be bfs");

        String json = pair.forward.toJson();
        assertTrue(json.contains("\"distance\":0"), "JSON must report distance 0");
        assertTrue(json.contains("\"solver\":\"bfs\""), "JSON must include solver bfs");
    }

    private static void testSingleRotation() {
        Triangulation base = new Triangulation(4, new Random(7));
        Triangulation target = base.clone2();

        boolean rotated = false;
        List<Integer> ids = target.getInternalNodeIds();
        for (int id : ids) {
            if (target.canRotateLeftById(id)) {
                target.rotateLeftById(id);
                rotated = true;
                break;
            }
            if (target.canRotateRightById(id)) {
                target.rotateRightById(id);
                rotated = true;
                break;
            }
        }
        assertTrue(rotated, "Expected at least one valid rotation for test case");

        TriangulationMetrics.Options opts = new TriangulationMetrics.Options();
        opts.timeLimitSeconds = 1.0;
        opts.visitedCap = 200_000;
        opts.queueCap = 200_000;

        TriangulationMetrics.ResultPair pair = TriangulationMetrics.runBidirectional(base, target, opts);
        assertTrue(pair.distanceAgrees, "Distances should agree in both directions");
        assertTrue(pair.forward.distance > 0, "Rotation distance must be positive");
        assertTrue("ok".equals(pair.forward.status), "Forward status should be ok");
        assertTrue("ok".equals(pair.reverse.status), "Reverse status should be ok");
        assertTrue("bfs".equals(pair.forward.solver), "Forward solver should be bfs");
        assertTrue("bfs".equals(pair.reverse.solver), "Reverse solver should be bfs");
    }

    private static void testCanonicalRoundTrip() {
        Triangulation base = new Triangulation(5, new Random(9));
        String canonical = base.canonicalTraversal();
        Triangulation rebuilt = Triangulation.fromCanonicalTraversal(canonical);

        TriangulationMetrics.Options opts = new TriangulationMetrics.Options();
        opts.timeLimitSeconds = 1.0;
        opts.visitedCap = 100_000;
        opts.queueCap = 100_000;

        TriangulationMetrics.ResultPair pair = TriangulationMetrics.runBidirectional(base, rebuilt, opts);
        assertTrue(pair.distanceAgrees, "Canonical rebuild should agree in both directions");
        assertTrue(pair.forward.distance == 0, "Canonical rebuild distance must be zero");
        assertTrue(canonical.equals(rebuilt.canonicalTraversal()), "Canonical traversal should round-trip");
        assertTrue("bfs".equals(pair.forward.solver), "Canonical solver should be bfs");
    }

    public static void main(String[] args) {
        testIdentity();
        testSingleRotation();
        testCanonicalRoundTrip();
        System.out.println("TriangulationMetricsTest OK");
    }
}
