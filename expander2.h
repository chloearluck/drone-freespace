#include <vector>
#include <map>
using namespace std;

class IloCplex;
class IloNumVarArray;

/* 
   For each vertex a involved in the LP, call

   addDisplacement(vert, (a - a_orig)/delta);

   where vert is the index of a and a_orig is its original position.

   For each pair of features A and B to be separated
   let p in A and q in B be the closest pair of points.  Declare

   Expander2::Pair pair((q-p).unit(), (q-p).length()/delta);

   Let u,v,w be orthonormal.

   For each a in A vertex, call

   pair.addConstraint(0, Constraint(vert, (a-p).dot(v)/delta, 
                                          (a-p).dot(w)/delta,
                                          (a-p).dot(u)/delta));

   where vert is the index of a.

   For each b in B vertex, call

   pair.addConstraint(1, Constraint(vert, (b-q).dot(v)/delta, 
                                          (b-q).dot(w)/delta,
                                          (b-q).dot(u)/delta));

   where vert is the index of b.

   Call expander2.addPair(pair).

   After adding all pairs, call expander2.expandV();

   The velocity of the vertex with index vert is

   expander2.getMotion(vert)

   Multiply by delta to get the actual velocity.
 */
class Expander2 {
public:
  class Point {
  public:
    double x[3];
    Point () { x[0] = x[1] = x[2] = 0; }
    Point (double a, double b, double c) {
      x[0] = a; x[1] = b; x[2] = c;
    }
    Point (double *p) {
      x[0] = p[0]; x[1] = p[1]; x[2] = p[2];
    }
    double dot (const Point &b) {
      return x[0] * b.x[0] + x[1] * b.x[1] + x[2] * b.x[2];
    }
  };

  class Constraint {
  public:
    int i;	// vertex index of a or b
    double v;	// (a-p)*v/delta or (b-q)*v/delta
    double w;	// (a-p)*w/delta or (b-q)*w/delta
    double r;	// (a-p)*u/delta or (b-q)*u/delta
    
    Constraint (int i, double v, double w, double r) 
      : i(i), v(v), w(w), r(r) {}
  };

  class Pair {
  public:
    Point u, v, w;
    double d;
    vector<Constraint> constraints[2]; // [0] a constraints, [1] b constraints

    // u is unit(q - p).  d is |q-p|/delta
    Pair (const Point &u, const Point &v, const Point &w, double d) : u(u), v(v), w(w), d(d) {}

    void addConstraint (int ab, const Constraint &constraint) {
      constraints[ab].push_back(constraint);
    }
  };

private:
  map<int, Point> disps;
  vector<Pair> pairs;
  map<int, Point> values;
  double uscaleBase;

  double expandV2 (double e, bool velocityObjective, double velocityBound, double uscale);

  bool checkPair (IloCplex &cplex, IloNumVarArray& cols,
                  map<int, int> &index,
                  int ipair, double t, double s);

public:
  Expander2 (double uscale=1.0) : uscaleBase(uscale) {}

  // vert is index of vertex and disp = (a_current - a_orig)/delta
  void addDisplacement (int vert, Point disp) { disps[vert] = disp; }

  void addPair (const Pair &pair) { pairs.push_back(pair); }

  bool expand ();

  // Expand by e*delta.
  // Motion is bounded by 1 (delta).
  bool expand (double e=1.0);

  // Expand by e*delta.
  bool expandV (double e=1.0, bool velocityObjective=true, double velocityBound=1.0);

  Point motion (int i) {
    map<int, Point>::iterator it = values.find(i);
    if (it == values.end())
      return Point(0, 0, 0);
    else
      return it->second;
  }
};
