import acm.program.*;
// import acm.graphics.*;
//import acm.util.*;
//import acm.gui.*;
//import java.util.Arrays;
import acm.program.GraphicsProgram;
import acm.gui.IntField;
import acm.graphics.GPoint;
import java.awt.*;
import java.util.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

/**
 * @Ge Xia, Copy Right 2017
 */

/**
 * TriangulationFlip.java <p>
 * 
 */
public class TriangulationFlip extends GraphicsProgram {
    // Counters for agreement across runs
    private int compareTotal = 0;
    private int compareEqual = 0;
    private JLabel agreeLabel;
    // initial size of the window

    public static int 
    APPLICATION_WIDTH = 1000,
    APPLICATION_HEIGHT = 1000;

    // instance variables 
    private int n = 10; // number of nodes, default is 10

    private int[] currentN = new int[2]; // the current triangulation displayed
    private Triangulation[][] triangulation = new Triangulation[2][]; // keep all previous triangulations in an array
    private int[] triangulationN = new int[2]; // number of triangulations in the array

    private Triangulation[] solution = new Triangulation[2]; 

    private Random rand = new Random();

    private IntField numField; // number of vertices
    private JButton generateButton, solveButton; // generate triangulation

    private JButton[] forwardButton = new JButton[2]; // move to next triangulation
    private JButton[] backwardButton = new JButton[2]; // move to previous triangulation
    private JButton[] fixButton = new JButton[2]; // fix the graph

    private JCheckBox skipSimpleCases;
    private boolean skipSimple = false;

    private byte[] dist = new byte[(int) Math.pow(2,31)-2]; // distance from s, is really dist + 1
    //private byte[][] dist2 = new byte[3][(int) Math.pow(2,31)-2]; // distance from s, is really dist + 1
    private Queue<Triangulation> bfsQueue = new LinkedList<Triangulation>();

    private void recordAgreement(boolean equal, int d01, int d10) {
        compareTotal++;
        if (equal) compareEqual++;
        int pct = (compareTotal == 0) ? 0 : (int)Math.round(100.0 * compareEqual / compareTotal);
        agreeLabel.setText("Matches: " + compareEqual + "/" + compareTotal + " (" + pct + "%)");
        System.out.println("A→B=" + d01 + ", B→A=" + d10 + ", equal=" + equal);
    }

    /** the run method, draw the inital graphics */
    public void init() {
        drawGUI();
    }

    private void drawGUI() {
        // forward button
        backwardButton[0] = new JButton("<");
        backwardButton[0].addActionListener(this);
        add(backwardButton[0], NORTH);

        // backward button
        forwardButton[0] = new JButton(">");
        forwardButton[0].addActionListener(this);
        add(forwardButton[0], NORTH);

        // fix button
        fixButton[0] = new JButton("Fix");
        fixButton[0].addActionListener(this);
        add(fixButton[0], NORTH);

        // status label
        add(new JLabel("Size: "), NORTH);

        // number of vertices to generate
        numField = new IntField(10);
        numField.setColumns(10);
        add(numField, NORTH);

        // generate button
        generateButton = new JButton("Generate");
        generateButton.addActionListener(this);
        add(generateButton, NORTH);

        // skip simple cases
        skipSimpleCases = new JCheckBox("Skip Simple Cases");
        add(skipSimpleCases, NORTH);

        // solve button
        solveButton = new JButton("Solve");
        solveButton.addActionListener(this);
        add(solveButton, NORTH);

        // forward button
        backwardButton[1] = new JButton("<");
        backwardButton[1].addActionListener(this);
        add(backwardButton[1], NORTH);

        // backward button
        forwardButton[1] = new JButton(">");
        forwardButton[1].addActionListener(this);
        add(forwardButton[1], NORTH);

        // fix button
        fixButton[1] = new JButton("Fix");
        fixButton[1].addActionListener(this);
        add(fixButton[1], NORTH);

        // enable or disable forward and backward buttons
        updateForwardBackwardButtons();

        agreeLabel = new JLabel("Matches: 0/0");
        add(agreeLabel, NORTH);
    }

    public void actionPerformed (ActionEvent evt) {
        if (evt.getSource() == generateButton) {  

            generateTriangulation();

            // regen until not simple
            skipSimple = skipSimpleCases.isSelected();
            int regenN = 0;
            if (skipSimple) {
                while (isSimpleCase()) {
                    generateTriangulation();
                    regenN ++;
                }
            }

            System.out.println("regenN: " + regenN);

        } else if (evt.getSource() == forwardButton[0]) {

            currentN[0] ++; // select the next triangulation as the current triangulation

        } else if (evt.getSource() == backwardButton[0]) {      
            currentN[0] --; // select the previous triangulation as the current triangulation

        } else if (evt.getSource() == fixButton[0]) {      

            fixTriangulation(0); // remove the history
        } else if (evt.getSource() == forwardButton[1]) {

            currentN[1] ++; // select the next triangulation as the current triangulation

        } else if (evt.getSource() == backwardButton[1]) {      

            currentN[1] --; // select the previous triangulation as the current triangulation

        } else if (evt.getSource() == fixButton[1]) {      

            fixTriangulation(1); // remove the history

        } else if (evt.getSource() == solveButton) {      
            int d01 = solveDistance(0, 1);
            int d10 = solveDistance(1, 0);
            boolean equal = (d01 == d10) && (d01 >= 0);
            recordAgreement(equal, d01, d10);
        }

        redraw();

    } 

    // check if it is a simple case, i.e., if there is a common edge or an easy flip
    private boolean isSimpleCase() {
        if (triangulationN[0] <= 0 || triangulationN[1] <= 0) {
            return true;
        }

        // source
        Triangulation s = triangulation[0][currentN[0]-1];
        Triangulation d = triangulation[1][currentN[1]-1];   

        Edge[] line = s.getLines();
        int lineN = s.getLineN();
        // check every edge in s, return true if there is a common edge or if the new edge after a flip is in d
        for (int i = 0; i < lineN; i++) {
            if (d.containsEdge(line[i])) { // it is simple if the edge is already in d
                return true;
            }

            Triangulation s1 = s.clone2();
            Edge newLine = s1.flip(i);
            if (d.containsEdge(newLine)) { // it is simple if the edge after flip is in d
                return true;
            }
        }

        // check every edge in d, return true if the new edge after a flip is in s
        line = d.getLines();
        lineN = d.getLineN();
        for (int i = 0; i < lineN; i++) {
            Triangulation d1 = d.clone2();
            Edge newLine = d1.flip(i);
            if (s.containsEdge(newLine)) { // it is simple if the edge after flip is in s
                return true;
            }
        }

        return false;
    }

    /*
    // solve the flip distance using a BFS search, transform from index start to index end (0 or 1), to produce index solN (0 or 1)
    private void solve(int start, int end, int solN) {
    if (triangulationN[0] <= 0 || triangulationN[1] <= 0) {
    return;
    }

    //System.out.print('\u000C');
    Arrays.fill(dist, (byte)0);
    bfsQueue.clear();

    // source
    Triangulation s = triangulation[start][currentN[start]-1];
    Triangulation d = triangulation[end][currentN[end]-1];

    // set dist of s, 0 is infinity, 1 is really 0
    dist[s.getRepN()] = 1;

    // add s to queue
    bfsQueue.add(s);

    Triangulation t = bfsQueue.poll();

    long steps = 0;

    // process the queue in BFS
    while (t != null) {
    //System.out.print(dist[t.getRepN()] + " -> ");
    if (t.getRepN() == d.getRepN()) {
    solution[solN] = t;
    //System.out.println("flips: " + (dist[t.getRepN()]-1));
    //System.out.println("steps: " + steps);
    return;
    }
    //System.out.println(Integer.toBinaryString(t.getRepN()) );
    //System.out.println(dist[t.getRepN()]);
    Triangulation[] candidate = new Triangulation[n];
    int cnt = 0;

    Edge[] line = t.getLines();
    int lineN = t.getLineN();
    // try every possible flips, do not flip edges already in d
    // if the new edge after a flip is in d, then this is the only candidate to consider
    for (int i = 0; i < lineN; i++) {
    if (!d.containsEdge(line[i])) { // only flip if the edge is not already in d
    Triangulation t1 = t.clone2();
    Edge newLine = t1.flip(i);
    if (d.containsEdge(newLine)) { // if the edge after flip is in d, this is the step to take
    if (dist[t1.getRepN()] == 0) { // not seen before
    candidate[0] = t1;
    cnt = 1;
    } else { // if seen before stop this branch
    cnt = 0;
    }
    break;
    } else { // otherwise, add the new triangulation to candidate if it is not already seen before
    if (dist[t1.getRepN()] == 0) { // not seen before
    candidate[cnt] = t1;
    cnt++;
    }
    }
    }
    }

    // add candidates to the queue in order
    for (int i = 0; i < cnt; i++) {
    dist[candidate[i].getRepN()] = (byte) (dist[t.getRepN()] + 1);
    //System.out.print(dist[candidate[i].getRepN()] + " | ");
    bfsQueue.add(candidate[i]);
    }

    //System.out.println(cnt);
    steps += cnt;

    // get the next triangulation on the queue
    t = bfsQueue.poll();
    }

    }
     */
    // fix the current triangulation as the only triangulation
    private void fixTriangulation(int k) {
        // remove all flipped edges and color the newest edge black
        triangulation[k][currentN[k]-1].fix();

        // create a clone of the current triangulation as the only triangulation
        //triangulation[k][0] = triangulation[k][currentN[k]-1].clone2();
        triangulation[k][0] = triangulation[k][currentN[k]-1];
        currentN[k] = 1;
        triangulationN[k] = 1;

        // redraw the triangulation
        redraw();

    }

    /** create a bouncy ball when mouse is pressed */
    public void mousePressed(GPoint point) {
        for (int k = 0; k < 2; k++) {
            if (currentN[k] > 0) {
                int i = triangulation[k][currentN[k]-1].findEdge(point);
                if (i >= 0) {
                    triangulation[k][currentN[k]] = triangulation[k][currentN[k]-1].clone2();
                    currentN[k]++;
                    triangulationN[k] = currentN[k];
                    triangulation[k][triangulationN[k]-1].flip(i);
                    redraw();
                } 
            }
        }

    }

    private void redraw() {
        removeAll();

        // Position the two binary trees side-by-side
        add(triangulation[0][currentN[0]-1], APPLICATION_WIDTH/4, APPLICATION_HEIGHT/4);
        add(triangulation[1][currentN[1]-1], APPLICATION_WIDTH*3/4, APPLICATION_HEIGHT/4);

        updateForwardBackwardButtons();
    }


    // update the Forward and Backward Buttons
    public void updateForwardBackwardButtons() {
        for (int k = 0; k < 2; k++) {
            forwardButton[k].setEnabled(currentN[k] != triangulationN[k]);

            backwardButton[k].setEnabled(currentN[k] > 1);
        }
    }

    // initially draw the walls and a label
    private void generateTriangulation() {

        n = numField.getValue();
        //System.out.println("generating two triangulations of " + n + " points");

        for (int k = 0; k < 2; k++) {
            triangulation[k] = new Triangulation[n*n];
            triangulationN[k] = 0;
            currentN[k] = triangulationN[k];

            triangulation[k][currentN[k]] = new Triangulation(n, rand);
            triangulation[k][currentN[k]].randTraingulate();
            triangulationN[k] = currentN[k] + 1;
            currentN[k] = triangulationN[k];
        }

        redraw();
    }

    /*
    // Return the minimal number of flips from startIndex to endIndex, or -1 if not found.
    private int solveDistance(int start, int end) {
    if (triangulationN[0] <= 0 || triangulationN[1] <= 0) return -1;

    Arrays.fill(dist, (byte)0);
    bfsQueue.clear();

    Triangulation s = triangulation[start][currentN[start]-1];
    Triangulation d = triangulation[end][currentN[end]-1];

    dist[s.getRepN()] = 1;      // store distance+1
    bfsQueue.add(s);

    while (!bfsQueue.isEmpty()) {
    Triangulation t = bfsQueue.poll();
    int tKey = t.getRepN();

    if (tKey == d.getRepN()) {
    return (dist[tKey] & 0xFF) - 1; // convert byte to int, subtract the +1
    }

    // generate candidates (same logic you already use)
    Triangulation[] candidate = new Triangulation[n];
    int cnt = 0;

    Edge[] line = t.getLines();
    int lineN = t.getLineN();
    for (int i = 0; i < lineN; i++) {
    if (!d.containsEdge(line[i])) { // only flip if not already in target
    Triangulation t1 = t.clone2();
    Edge newLine = t1.flip(i);
    int key = t1.getRepN();
    if (d.containsEdge(newLine)) {
    if (dist[key] == 0) {
    candidate[0] = t1;
    cnt = 1;
    } else {
    cnt = 0;
    }
    break;
    } else {
    if (dist[key] == 0) {
    candidate[cnt++] = t1;
    }
    }
    }
    }

    for (int i = 0; i < cnt; i++) {
    int key = candidate[i].getRepN();
    dist[key] = (byte)((dist[tKey] & 0xFF) + 1);
    bfsQueue.add(candidate[i]);
    }
    }

    return -1;
    }
     */

    // Compute the minimum flip (rotation) distance between triangulation[start] and triangulation[end]
private int solveDistance(int start, int end) {
    if (triangulationN[start] <= 0 || triangulationN[end] <= 0) return -1;

    Triangulation source = triangulation[start][currentN[start] - 1];
    Triangulation target = triangulation[end][currentN[end] - 1];

    final String targetRep = target.getRep();   // FIXED target key

    if (source.getRep().equals(targetRep)) return 0;

    Queue<Triangulation> queue = new LinkedList<>();
    Map<String, Integer> distMap = new HashMap<>();

    queue.add(source);
    distMap.put(source.getRep(), 0);

    while (!queue.isEmpty()) {
        Triangulation cur = queue.poll();
        int curDist = distMap.get(cur.getRep());

        for (int id : cur.getInternalNodeIds()) {
            if (cur.canRotateRightById(id)) {
                Triangulation nxt = cur.clone2();
                nxt.rotateRightById(id);
                String k = nxt.getRep();
                if (!distMap.containsKey(k)) {
                    if (k.equals(targetRep)) return curDist + 1;
                    distMap.put(k, curDist + 1);
                    queue.add(nxt);
                }
            }
            if (cur.canRotateLeftById(id)) {
                Triangulation nxt = cur.clone2();
                nxt.rotateLeftById(id);
                String k = nxt.getRep();
                if (!distMap.containsKey(k)) {
                    if (k.equals(targetRep)) return curDist + 1;
                    distMap.put(k, curDist + 1);
                    queue.add(nxt);
                }
            }
        }
    }
    return -1;
}


    public static void main(String[] args) {
        // This class is mandatory to be executed by the JVM.
        // You don't need to do anything in here, since you're subclassing ConsoleProgram,
        // which invokes the run() method.
        new TriangulationFlip().start(args);
    }

}