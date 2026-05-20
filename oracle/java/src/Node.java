/**
 * @Ge Xia, Copy Right 2017
 */

import acm.graphics.*;

public class Node extends GCompound {
    private int index;

    public Node(int index, double w, double h) {
        this.index = index;
        GOval oval = new GOval(w, h);
        oval.setFilled(true);
        add(oval, -w/2, -h/2);
    }

    public Node(int index, double x, double y, double w, double h) {
        this(index, w, h);
        setLocation(x, y);
    }
    
    public int getIndex() {
      return index;
    }
}