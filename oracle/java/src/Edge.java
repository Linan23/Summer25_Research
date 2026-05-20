/**
 * @Ge Xia, Copy Right 2017
 */

import acm.graphics.*;
import java.awt.*;


/*
public class Edge extends GCompound {
  private int startIndex, endIndex;
  private GLine line;

  // startIndex is always less than endIndex
  public Edge(int startIndex, int endIndex, double x1, double y1, 
  double x2, double y2, Triangulation triangulaton) {
    this.startIndex = startIndex;
    this.endIndex = endIndex;
    line = new GLine(x1, y1, x2, y2);
    add(line);
    //triangulaton.add(line);
  }

  public GPoint getStartPoint() {
    return line.getStartPoint();
  }

  public GPoint getEndPoint() {
    return line.getEndPoint();
  }

  public boolean contains2(GPoint p) {
    GPoint p1 = new GPoint(p.getX()+1, p.getY());
    GPoint p2 = new GPoint(p.getX()-1, p.getY());
    GPoint p3 = new GPoint(p.getX(), p.getY()+1);
    GPoint p4 = new GPoint(p.getX(), p.getY()-1);
    
    return line.contains(p) || line.contains(p1) || line.contains(p2) || 
    line.contains(p3) || line.contains(p4);
  }

  public void setColor(Color color) {
    line.setColor(color);
  }
  
  public int getStartIndex() {
    return startIndex;
  }
  
  public int getEndIndex() {
    return endIndex;
  }
  
  public Color getColor() {
    return line.getColor();
  }
  
  public boolean match(int i, int j) {
    return (startIndex == i && endIndex == j) || (startIndex == j && endIndex == i);
  }
}
*/

public class Edge {
    public int id1, id2;

    public Edge(int i, int j) {
        this.id1 = i;
        this.id2 = j;
    }

    // equality test
    public boolean match(int a, int b) {
        return (id1 == a && id2 == b) || (id1 == b && id2 == a);
    }
}

