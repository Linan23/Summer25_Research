import java.util.*;
import acm.graphics.*;

/*
public class Triangulation {
private TreeNode root;
private int n;
private Random rand;
private String rep;
private int repN;

// Tree node structure (internal node has left/right children)
private static class TreeNode {
TreeNode left, right, parent;
}

public Triangulation(int n, Random rand) {
this.n = n;
this.rand = rand;
this.root = buildRandomTree(n, rand);
computeRep();
}

// Random full binary tree generator (n internal nodes)
private TreeNode buildRandomTree(int n, Random rand) {
if (n <= 0) return null;
TreeNode node = new TreeNode();
int leftSize = rand.nextInt(n);  // random split
node.left = buildRandomTree(leftSize, rand);
if (node.left != null) node.left.parent = node;
node.right = buildRandomTree(n - 1 - leftSize, rand);
if (node.right != null) node.right.parent = node;
return node;
}

public Triangulation clone2() {
Triangulation t = new Triangulation(0, rand);
t.root = cloneTree(this.root);
t.n = this.n;
t.computeRep();
return t;
}

private TreeNode cloneTree(TreeNode node) {
if (node == null) return null;
TreeNode c = new TreeNode();
c.left = cloneTree(node.left);
if (c.left != null) c.left.parent = c;
c.right = cloneTree(node.right);
if (c.right != null) c.right.parent = c;
return c;
}

// Perform rotation (flip) on a given internal edge
public void flipRandomEdge() {
List<TreeNode> internal = collectInternal(root);
if (internal.isEmpty()) return;
TreeNode pivot = internal.get(rand.nextInt(internal.size()));
rotate(pivot);
computeRep();
}

private void rotate(TreeNode parent) {
if (parent == null) return;
// randomly choose left or right rotation
boolean left = rand.nextBoolean();
if (left && parent.left != null)
rotateRight(parent);
else if (!left && parent.right != null)
rotateLeft(parent);
}

private void rotateLeft(TreeNode x) {
TreeNode y = x.right;
if (y == null) return;
x.right = y.left;
if (y.left != null) y.left.parent = x;
y.parent = x.parent;
if (x.parent == null) root = y;
else if (x.parent.left == x) x.parent.left = y;
else x.parent.right = y;
y.left = x;
x.parent = y;
}

private void rotateRight(TreeNode x) {
TreeNode y = x.left;
if (y == null) return;
x.left = y.right;
if (y.right != null) y.right.parent = x;
y.parent = x.parent;
if (x.parent == null) root = y;
else if (x.parent.left == x) x.parent.left = y;
else x.parent.right = y;
y.right = x;
x.parent = y;
}

private List<TreeNode> collectInternal(TreeNode node) {
List<TreeNode> res = new ArrayList<>();
if (node == null) return res;
if (node.left != null || node.right != null) res.add(node);
res.addAll(collectInternal(node.left));
res.addAll(collectInternal(node.right));
return res;
}

// Compute bracket/Dyck-word rep (same as original Triangulation)
private void computeRep() {
StringBuilder sb = new StringBuilder();
encode(root, sb);
rep = sb.toString().replace("a", "");
repN = Integer.parseInt(rep.replace("a", ""), 2);
}

private void encode(TreeNode node, StringBuilder sb) {
if (node == null) return;
sb.append("0"); // open
encode(node.left, sb);
sb.append("a");
encode(node.right, sb);
sb.append("1"); // close
}

public void randTraingulate() {
// no-op, since already randomized tree structure
}

public void fix() {}  // keep placeholder

public String getRep() { return rep; }
public int getRepN() { return repN; }

public boolean IsEquivalent(Triangulation t) {
return this.rep.equals(t.rep);
}

public static class Edge {
int id1, id2;
public Edge(int i, int j) { id1 = i; id2 = j; }
}

// pretend that each internal node represents an "edge"
public Edge[] getLines() {
List<Edge> list = new ArrayList<>();
List<TreeNode> internal = collectInternal(root);
int counter = 0;
for (TreeNode t : internal)
list.add(new Edge(counter++, counter + 1));
return list.toArray(new Edge[0]);
}

public int getLineN() {
return collectInternal(root).size();
}

public Edge flip(int index) {
List<TreeNode> internal = collectInternal(root);
if (index < 0 || index >= internal.size()) return new Edge(-1, -1);
TreeNode target = internal.get(index);
if (rand.nextBoolean()) rotateLeft(target);
else rotateRight(target);
computeRep();
return new Edge(index, index + 1);
}

public boolean containsEdge(Edge e) {
return false;
}
}
 */

public class Triangulation extends GCompound {
    private TreeNode root;
    private int n;
    private Random rand;
    private String rep;
    private int repN;
    private int nextId = 1; // Unique ID counter for nodes (monotone increasing per tree)
    
    // drawing constants
    private static final double NODE_RADIUS = 12;
    private static final double X_SPACING = 40;
    private static final double Y_SPACING = 50;

    // internal binary tree node
    private static class TreeNode {
        TreeNode left, right, parent;
        int id;
        double x, y; // coordinates for drawing
    }

    // constructor
    public Triangulation(int n, Random rand) {
        this.n = n;
        this.rand = rand;
        this.nextId = 1;                 // reset counter for a fresh tree
        this.root = buildRandomTree(n);
        assertUniqueIds();
        computeRep();
        layoutAndDraw();
    }

    private TreeNode buildRandomTree(int n) {
    if (n <= 0) return null;
    TreeNode node = new TreeNode();
    node.id = nextId++;                 // ← UNIQUE, STABLE ID

    // keep the random shape; only IDs must be deterministic & unique
    int leftSize = (n == 1) ? 0 : rand.nextInt(n);  // 0..n-1
    node.left = buildRandomTree(leftSize);
    if (node.left != null) node.left.parent = node;
    node.right = buildRandomTree(n - 1 - leftSize);
    if (node.right != null) node.right.parent = node;
    return node;
    }
    
    /*
    public Triangulation clone2() {
        Triangulation t = new Triangulation(0, rand);
        t.root = cloneTree(this.root);
        t.n = this.n;
        t.computeRep();
        t.layoutAndDraw();
        return t;
    }
    */
       
    public Triangulation clone2() {
    Triangulation t = new Triangulation(0, rand);
    t.root = cloneTree(this.root);
    t.n = this.n;
    t.nextId = Math.max(this.nextId, 1 + maxId(t.root));
    t.assertUniqueIds();
    t.computeRep();
    return t;
}

private int maxId(TreeNode node) {
    if (node == null) return 0;
    return Math.max(node.id, Math.max(maxId(node.left), maxId(node.right)));
}

    private TreeNode cloneTree(TreeNode node) {
        if (node == null) return null;
        TreeNode c = new TreeNode();
        c.id = node.id;
        c.left = cloneTree(node.left);
        if (c.left != null) c.left.parent = c;
        c.right = cloneTree(node.right);
        if (c.right != null) c.right.parent = c;
        return c;
    }

    private void assertUniqueIds() {
    java.util.HashSet<Integer> seen = new java.util.HashSet<>();
    for (int id : getInternalNodeIds()) {
        if (!seen.add(id)) {
            throw new IllegalStateException("Duplicate node id: " + id);
        }
    }
}

    private static String joinIds(List<Integer> values) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < values.size(); ++i) {
            if (i > 0) sb.append(',');
            sb.append(values.get(i));
        }
        return sb.toString();
    }

    private void collectPreorder(TreeNode node, List<Integer> out) {
        if (node == null) return;
        out.add(node.id);
        collectPreorder(node.left, out);
        collectPreorder(node.right, out);
    }

    private void collectInorder(TreeNode node, List<Integer> out) {
        if (node == null) return;
        collectInorder(node.left, out);
        out.add(node.id);
        collectInorder(node.right, out);
    }

    public String canonicalTraversal() {
        List<Integer> pre = new ArrayList<>();
        List<Integer> in  = new ArrayList<>();
        collectPreorder(root, pre);
        collectInorder(root, in);
        return "P:" + joinIds(pre) + ";I:" + joinIds(in);
    }

    public static Triangulation fromTraversals(List<Integer> preorder, List<Integer> inorder) {
        if (preorder == null || inorder == null)
            throw new IllegalArgumentException("Traversals must not be null");
        if (preorder.size() != inorder.size())
            throw new IllegalArgumentException("Preorder and inorder must have equal length");

        Triangulation t = new Triangulation(0, new Random(0));
        t.n = inorder.size();

        Map<Integer, Integer> index = new java.util.HashMap<>();
        for (int i = 0; i < inorder.size(); ++i) {
            index.put(inorder.get(i), i);
        }

        t.root = t.buildFromTraversals(preorder, 0, preorder.size(), inorder, 0, inorder.size(), index);
        t.nextId = Math.max(t.nextId, 1 + t.maxId(t.root));
        t.assertUniqueIds();
        t.computeRep();
        return t;
    }

    public static Triangulation fromCanonicalTraversal(String canonical) {
        if (canonical == null)
            throw new IllegalArgumentException("Canonical traversal string must not be null");
        if (!canonical.startsWith("P:"))
            throw new IllegalArgumentException("Canonical traversal must start with \"P:\"");
        int split = canonical.indexOf(";I:");
        if (split < 0)
            throw new IllegalArgumentException("Canonical traversal missing \";I:\" delimiter");

        String prePart = canonical.substring(2, split);
        String inPart = canonical.substring(split + 3);
        List<Integer> preorder = parseSequence(prePart);
        List<Integer> inorder = parseSequence(inPart);
        return fromTraversals(preorder, inorder);
    }

    private static List<Integer> parseSequence(String part) {
        List<Integer> out = new ArrayList<>();
        if (part == null || part.isEmpty()) return out;
        String[] tokens = part.split(",");
        for (String token : tokens) {
            if (token == null || token.isEmpty()) continue;
            out.add(Integer.parseInt(token.trim()));
        }
        return out;
    }

    private TreeNode buildFromTraversals(List<Integer> preorder, int ps, int pe,
                                         List<Integer> inorder, int is, int ie,
                                         Map<Integer, Integer> index) {
        if (ps >= pe || is >= ie) return null;
        int rootVal = preorder.get(ps);
        Integer pos = index.get(rootVal);
        if (pos == null)
            throw new IllegalArgumentException("Value " + rootVal + " missing from inorder traversal");
        if (pos < is || pos >= ie)
            throw new IllegalArgumentException("Invalid inorder position for value " + rootVal);

        TreeNode node = new TreeNode();
        node.id = rootVal;

        int leftSize = pos - is;
        node.left = buildFromTraversals(preorder, ps + 1, ps + 1 + leftSize,
                                        inorder, is, pos, index);
        if (node.left != null) node.left.parent = node;

        node.right = buildFromTraversals(preorder, ps + 1 + leftSize, pe,
                                         inorder, pos + 1, ie, index);
        if (node.right != null) node.right.parent = node;

        return node;
    }
    // ----------------------------------------------------------------------
    // Rotation logic
    // ----------------------------------------------------------------------
    private void rotateLeft(TreeNode x) {
        TreeNode y = x.right;
        if (y == null) return;
        x.right = y.left;
        if (y.left != null) y.left.parent = x;
        y.parent = x.parent;
        if (x.parent == null) root = y;
        else if (x.parent.left == x) x.parent.left = y;
        else x.parent.right = y;
        y.left = x;
        x.parent = y;
    }
    
    private void rotateRight(TreeNode x) {
        TreeNode y = x.left;
        if (y == null) return;
        x.left = y.right;
        if (y.right != null) y.right.parent = x;
        y.parent = x.parent;
        if (x.parent == null) root = y;
        else if (x.parent.left == x) x.parent.left = y;
        else x.parent.right = y;
        y.right = x;
        x.parent = y;
    }
    
    
    private TreeNode findById(TreeNode node, int id) {
    if (node == null) return null;
    if (node.id == id) return node;
    TreeNode L = findById(node.left, id);
    return (L != null) ? L : findById(node.right, id);
}

// expose safe checks
public boolean canRotateLeftById(int id) {
    TreeNode x = findById(root, id);
    return x != null && x.right != null;
}

public boolean canRotateRightById(int id) {
    TreeNode x = findById(root, id);
    return x != null && x.left != null;
}

// pure rotations for solver (no layout/draw)
public void rotateLeftById(int id) {
    TreeNode x = findById(root, id);
    if (x != null && x.right != null) {
        rotateLeft(x);
        computeRep();
    }
}

public void rotateRightById(int id) {
    TreeNode x = findById(root, id);
    if (x != null && x.left != null) {
        rotateRight(x);
        computeRep();
    }
}
    

    private List<TreeNode> collectInternal(TreeNode node) {
        List<TreeNode> res = new ArrayList<>();
        if (node == null) return res;
        if (node.left != null || node.right != null) res.add(node);
        res.addAll(collectInternal(node.left));
        res.addAll(collectInternal(node.right));
        return res;
    }
    
    public List<Integer> getInternalNodeIds() {
    List<Integer> ids = new ArrayList<>();
    for (TreeNode t : collectInternal(root)) ids.add(t.id);
    java.util.Collections.sort(ids);   // deterministic expansion order
    return ids;
}
    // ----------------------------------------------------------------------
    // Compatibility layer (same API as before)
    // ----------------------------------------------------------------------
    public Edge[] getLines() {
        List<Edge> list = new ArrayList<>();
        List<TreeNode> internal = collectInternal(root);
        int counter = 0;
        for (TreeNode t : internal)
            list.add(new Edge(counter, counter + 1));
        return list.toArray(new Edge[0]);
    }

    public int getLineN() {
        return collectInternal(root).size();
    }

    // UI helper: flip left/right explicitly; keeps drawing for the UI only
public Edge flip(int index, boolean toLeft) {
    List<TreeNode> internal = collectInternal(root);
    if (index < 0 || index >= internal.size()) return new Edge(-1, -1);
    TreeNode target = internal.get(index);

    if (toLeft && target.right != null) rotateLeft(target);
    else if (!toLeft && target.left != null) rotateRight(target);

    computeRep();
    layoutAndDraw(); // drawing is UI-only
    return new Edge(index, index + 1);
}

public Edge flip(int index) {
    List<TreeNode> internal = collectInternal(root);
    if (index < 0 || index >= internal.size()) return new Edge(-1, -1);
    TreeNode target = internal.get(index);
    boolean toLeft = (target.right != null); // prefer left-rotation when possible
    return flip(index, toLeft);
}

    public boolean containsEdge(Edge e) { return false; }

    public int findEdge(GPoint point) {
        if (getLineN() == 0) return -1;
        return rand.nextInt(getLineN());
    }

    public void randTraingulate() {}

    public void fix() {}

    // ----------------------------------------------------------------------
    // Representation (binary encoding)
    // ----------------------------------------------------------------------
    private void computeRep() {
        StringBuilder sb = new StringBuilder();
        encodeCanonical(root, sb);
        this.rep = sb.toString(); 
    }

    private void encodeCanonical(TreeNode node, StringBuilder sb) {
    if (node == null) {
        sb.append('#');           // explicit null marker
        return;
    }
    sb.append('(');
    encodeCanonical(node.left, sb);
    sb.append(',');
    encodeCanonical(node.right, sb);
    sb.append(')');
}

    public String getRep() { return rep; }

    public int getRepN() { return repN; }

    public boolean IsEquivalent(Triangulation t) { return this.rep.equals(t.rep); }

    // ----------------------------------------------------------------------
    // Visual layout and drawing
    // ----------------------------------------------------------------------
    private void layoutAndDraw() {
        removeAll();  // clear GCompound contents
        double widthScale = Math.max(120, 600.0 / (n + 2));  // adaptive horizontal spread
        computeNodePositions(root, 0, 0, widthScale, Y_SPACING);
        drawTree(root);
    }

    private int depth(TreeNode node) {
        if (node == null) return 0;
        return 1 + Math.max(depth(node.left), depth(node.right));
    }

    private void computeNodePositions(TreeNode node, double x, double y, double offset, double yStep) {
        if (node == null) return;
        node.x = x;
        node.y = y;
        if (node.left != null)
            computeNodePositions(node.left, x - offset / 2, y + yStep, offset / 2, yStep);
        if (node.right != null)
            computeNodePositions(node.right, x + offset / 2, y + yStep, offset / 2, yStep);
    }

    private void drawTree(TreeNode node) {
        if (node == null) return;

        // draw edges first
        if (node.left != null) {
            GLine l = new GLine(node.x, node.y, node.left.x, node.left.y);
            add(l);
        }
        if (node.right != null) {
            GLine l = new GLine(node.x, node.y, node.right.x, node.right.y);
            add(l);
        }

        // draw node circle
        GOval circle = new GOval(node.x - NODE_RADIUS, node.y - NODE_RADIUS, NODE_RADIUS * 2, NODE_RADIUS * 2);
        circle.setFilled(true);
        add(circle);

        drawTree(node.left);
        drawTree(node.right);
    }
    
    public void recomputeRep() {
    computeRep();
}
}
