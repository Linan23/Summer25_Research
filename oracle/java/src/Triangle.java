/**
 * @Ge Xia, Copy Right 2017
 */

import acm.graphics.*;

public class Triangle {
  private int v1, v2, v3;


  // counter cc-wise, v1, e1, v2, e2, v3, e3
  public Triangle(int v1, int v2, int v3) {
    this.v1 = v1;
    this.v2 = v2;
    this.v3 = v3;
  }

  public boolean contains(Edge e) {
    return matchE1(e) || matchE2(e) || matchE3(e);
  }

  private boolean matchE1(Edge e) {
    return e.match(v1, v2);
  }

  private boolean matchE2(Edge e) {
    return e.match(v2, v3);
  }

  private boolean matchE3(Edge e) {
    return e.match(v3, v1);
  }

  // vertex opposing e
  public int oppose(Edge e) {
    if (matchE1(e)) return v3;
    else if (matchE2(e)) return v1;
    else return v2;
  }

  // start of edge cc-wise
  public int start(Edge e) {
    if (matchE1(e)) return v1;
    else if (matchE2(e)) return v2;
    else return v3;
  }

  // end of edge cc-wise
  public int end(Edge e) {
    if (matchE1(e)) return v2;
    else if (matchE2(e)) return v3;
    else return v1;
  }
  
  public int getV1() {
    return v1;
  }
  
  public int getV2() {
    return v2;
  }
  
  public int getV3() {
    return v3;
  }
}