#include "expander2.h"
#include <ilcplex/ilocplex.h>
#include <string>

void Expander2::Pair::reorient () {
  int abR, iR, nR = 0;

  for (int ab = 0; ab <= 1; ab++)
    for (int i = 0; i < constraints[ab].size(); i++)
      if (constraints[ab][i].r != 0) {
	abR = ab;
	iR = i;
	nR++;
      }

  if (nR != 1)
    return;

  int ab1, i1, ab2, i2[2];

  if (constraints[0].size() == 2) {
    assert(constraints[1].size() == 2);
    ab1 = abR;
    ab2 = !abR;
    i1 = !iR;
    i2[0] = 0;
    i2[1] = 1;
  }
  else {
    ab1 = !abR;
    ab2 = abR;
    assert(constraints[ab1].size() == 1);
    assert(constraints[ab2].size() == 3);
    i1 = 0;
    i2[0] = (iR+1)%3;
    i2[1] = (iR+2)%3;
  }
  
  int i2i = (fabs(constraints[ab2][i2[0]].v) < fabs(constraints[ab2][i2[1]].v));
  double x = constraints[ab2][i2[i2i]].v;
  double y = constraints[ab2][i2[i2i]].w;
  double s = sqrt(x * x + y * y);
  x /= s;
  y /= s;
  Point v_ = v * x + w * y;
  Point w_ = v * -y + w * x;
  v = v_;
  w = w_;

  double x1 = constraints[ab1][i1].v * x + constraints[ab1][i1].w * y;
  double y1 = constraints[ab1][i1].v * -y + constraints[ab1][i1].w * x;
  constraints[ab1][i1].v = x1;
  constraints[ab1][i1].w = y1;

  constraints[ab2][i2[i2i]].v = s;
  constraints[ab2][i2[i2i]].w = 0;

  x = constraints[ab2][i2[!i2i]].v;
  y = constraints[ab2][i2[!i2i]].w;
  s = sqrt(x * x + y * y);
  constraints[ab2][i2[!i2i]].v = -s;
  constraints[ab2][i2[!i2i]].w = 0;
}

class Ineq {
public:
  bool gt;
  double b;
  vector<int> cols;
  vector<double> coefs;

  Ineq (bool gt, double b) : gt(gt), b(b) {}

  void add (int col, double coef) { 
    if (coef != 0) {
      cols.push_back(col);
      coefs.push_back(coef);
    }
  }

  double value (IloCplex &cplex, IloNumVarArray& ccols) {
    double v = -b;
    for (int i = 0; i < cols.size(); i++) {
      // cerr << "extracting " << cols[i] << " " << ccols[cols[i]].getName() << endl;
      v += cplex.getValue(ccols[cols[i]]) * coefs[i];
    }
    // cerr << endl;
    return gt ? -v : v;
  }
};

bool Expander2::expand () {
  vector<Vertex*> verts;
  map<Vertex*, int> index;

  int npairs = pairs.size();

  for (int ipair = 0; ipair < npairs; ipair++)
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Vertex* vert = constraints[ivert].i;
        if (index.find(vert) == index.end()) {
          index[vert] = verts.size();
          verts.push_back(vert);
        }
      }
    }

  int nverts = verts.size();

  IloEnv env;
  IloNumVarArray cols(env);

  // a' or b' is ((xp - xm), (yp - ym), (zp - zm))
  string xyzpm[] = { "xp", "xm", "yp", "ym", "zp", "zm" };
  for (int i = 0; i < nverts; i++) {
    string base = "v";
    base += std::to_string(i);
    for (int j = 0; j < 6; j++) {
      cols.add(IloNumVar(env));
      // cerr << "cols.size() " << cols.length << " = " << 6*i+j << endl;
      cols[6*i+j].setName((base + xyzpm[j]).c_str());
    }
  }

  // p' is (x,y,z).  u' = (vp - vm) v + (wp - wm) w
  string xyzvw[] = { "x", "y", "z", "vp", "vm", "wp", "wm" };
  for (int i = 0; i < npairs; i++) {
    string base = "p";
    base += std::to_string(i);
    for (int j = 0; j < 3; j++) {
      cols.add(IloNumVar(env, -IloInfinity, IloInfinity));
      cols[6*nverts + 7*i+j].setName((base + xyzvw[j]).c_str());
    }
    for (int j = 3; j < 7; j++) {
      cols.add(IloNumVar(env, 0, 1));
      cols[6*nverts + 7*i+j].setName((base + xyzvw[j]).c_str());
    }
  }

  IloObjective obj = IloMinimize(env);

  for (int i = 0; i < nverts; i++)
    for (int j = 0; j < 6; j++)
      obj.setLinearCoef(cols[6*i+j], 1.0);

  bool uObj = false;
  if (uObj) {
    for (int i = 0; i < npairs; i++)
      for (int j = 3; j < 7; j++)
        obj.setLinearCoef(cols[6*nverts+7*i+j], 1.0);
  }

  IloRangeArray rows(env);

  int nrows = 0;
  for (int ipair = 0; ipair < npairs; ipair++) {
    string base = "p";
    base += std::to_string(ipair);
    Pair &pair = pairs[ipair];
    int pCol = 6 * nverts + 7 * ipair;
    
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pair.constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Constraint con = constraints[ivert];
        int vCol = 6 * index[con.i];
        if (ifeat == 0) {
          rows.add(IloRange(env, -IloInfinity, -con.r));
          rows[nrows].setName((base + "a" + std::to_string(index[con.i])).c_str());
        }
        else {
          rows.add(IloRange(env, 1 - con.r, IloInfinity));
          rows[nrows].setName((base + "b" + std::to_string(index[con.i])).c_str());
        }
        
        IloRange row = rows[nrows];

        // (a' - p') * u + x (a - p) * v + y (a - p) * w
        for (int xyz = 0; xyz < 3; xyz++) {
          for (int pm = 0; pm < 2; pm++)
            row.setLinearCoef(cols[vCol + 2*xyz + pm],
			      pair.u.x[xyz] * (1-2*pm));
          row.setLinearCoef(cols[pCol + xyz], -pair.u.x[xyz]);
        }

	// double d = 1.01e-6;
	double d = 1;

        row.setLinearCoef(cols[pCol + 3],  con.v / d);
        row.setLinearCoef(cols[pCol + 4], -con.v / d);
        row.setLinearCoef(cols[pCol + 5],  con.w / d);
        row.setLinearCoef(cols[pCol + 6], -con.w / d);

        nrows++;
      }
    }
  }

  IloModel model(env);

  model.add(obj);
  model.add(rows);

  IloCplex cplex(model);

  // cplex.exportModel("expand.lp");

  cplex.setOut(env.getNullStream());
  if ( !cplex.solve() ) {
    env.error() << "Failed to optimize LP" << endl;
    exit(-1);
    return false;
  }
  
  bool verbose = false;

  IloNumArray vals(env);
  if (verbose) env.out() << "Solution status = " << cplex.getStatus() << endl;
  if (verbose) env.out() << "Solution value  = " << cplex.getObjValue() << endl;
  // cplex.getValues(vals, cols);
  //env.out() << "Values        = " << vals << endl;

  for (int i = 0; i < nverts; i++) {
    double x[3] = { 0, 0, 0 };
    for (int xyz = 0; xyz < 3; xyz++)
      for (int pm = 0; pm < 2; pm++)
        // x[xyz] += vals[6*i + 2*xyz + pm] * (1-2*pm);
	x[xyz] += cplex.getValue(cols[6*i + 2*xyz + pm]) * (1-2*pm);
    values[verts[i]] = Point(x);

    if (verbose)
      cerr << verts[i] << " "
	   << values[verts[i]].x[0] << " "
	   << values[verts[i]].x[1] << " "
	   << values[verts[i]].x[2] << endl;
  }

  env.end();
  
  return true;
}

bool Expander2::expand (double e) {
  vector<Vertex*> verts;
  map<Vertex*, int> index;

  int npairs = pairs.size();

  for (int ipair = 0; ipair < npairs; ipair++)
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Vertex* vert = constraints[ivert].i;
        if (index.find(vert) == index.end()) {
          index[vert] = verts.size();
          verts.push_back(vert);
        }
      }
    }

  int nverts = verts.size();

  IloEnv env;
  IloNumVarArray cols(env);
  IloNumArray coefs(env);

  // a' or b' is ((xp - xm), (yp - ym), (zp - zm))
  string xyzpm[] = { "xp", "xm", "yp", "ym", "zp", "zm" };
  for (int i = 0; i < nverts; i++) {
    string base = "v";
    base += std::to_string(i);
    for (int j = 0; j < 6; j++) {
      cols.add(IloNumVar(env, 0, 1));
      // cerr << "cols.size() " << cols.length << " = " << 6*i+j << endl;
      cols[6*i+j].setName((base + xyzpm[j]).c_str());
      coefs.add(IloNum(1));
    }
  }

  // p' is (x,y,z).  u' = (vp - vm) v + (wp - wm) w
  string pvw[] = { "pp", "pm", "vp", "vm", "wp", "wm" };
  for (int i = 0; i < npairs; i++) {
    string base = "p";
    base += std::to_string(i);
    for (int j = 0; j < 6; j++) {
      if (j < 3) // ??? BUG should be j < 2 ?
	cols.add(IloNumVar(env, 0, sqrt(3.0)));
      else
	cols.add(IloNumVar(env, 0, 1));
      cols[6*nverts + 6*i+j].setName((base + pvw[j]).c_str());
      coefs.add(IloNum(1));
    }
  }

  int dCol = 6*nverts + 6*npairs;
  cols.add(IloNumVar(env, 0, e));
  cols[dCol].setName("D");
  coefs.add(IloNum(-6 * nverts));

  IloObjective obj = IloMinimize(env);

  if (false) {
    cerr << "setting objective coefficients of vertices" << endl;
    for (int i = 0; i < nverts; i++)
      for (int j = 0; j < 6; j++)
	obj.setLinearCoef(cols[6*i+j], 1.0);
    cerr << "done" << endl;
    
    cerr << "setting objective coefficients of pairs" << endl;
    bool uObj = true;
    if (uObj) {
      for (int i = 0; i < npairs; i++)
	for (int j = 0; j < 6; j++)
	  obj.setLinearCoef(cols[6*nverts+6*i+j], 1.0);
    }
    cerr << "done" << endl;
    
    obj.setLinearCoef(cols[dCol], -6 * nverts);
  }
  else {
    //cerr << "setting objective" << endl;
    obj.setLinearCoefs(cols, coefs);
    //cerr << "done" << endl;
  }
    
  IloRangeArray rows(env);

  int nrows = 0;
  for (int ipair = 0; ipair < npairs; ipair++) {
    string base = "p";
    base += std::to_string(ipair);
    Pair &pair = pairs[ipair];
    int pCol = 6 * nverts + 6 * ipair;
    
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pair.constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Constraint con = constraints[ivert];
        int vCol = 6 * index[con.i];
        if (ifeat == 0) {
          rows.add(IloRange(env, -IloInfinity, -con.r));
          rows[nrows].setName((base + "a" + std::to_string(index[con.i])).c_str());
        }
        else {
          rows.add(IloRange(env, - con.r, IloInfinity));
          rows[nrows].setName((base + "b" + std::to_string(index[con.i])).c_str());
        }
        
        IloRange row = rows[nrows];

        // (a' - p') * u + x (a - p) * v + y (a - p) * w
	for (int pm = 0; pm < 2; pm++) {
	  for (int xyz = 0; xyz < 3; xyz++)
            row.setLinearCoef(cols[vCol + 2*xyz + pm],
			      pair.u.x[xyz] * (1-2*pm));
	  row.setLinearCoef(cols[pCol + pm], -(1-2*pm));
        }

        row.setLinearCoef(cols[pCol + 2],  con.v);
        row.setLinearCoef(cols[pCol + 3], -con.v);
        row.setLinearCoef(cols[pCol + 4],  con.w);
        row.setLinearCoef(cols[pCol + 5], -con.w);

	if (ifeat == 1)
	  row.setLinearCoef(cols[dCol], -1);

        nrows++;
      }
    }
  }

  IloModel model(env);

  model.add(obj);
  model.add(rows);

  IloCplex cplex(model);

  static int iCall = 100;
  string base = "expand";
  // cplex.exportModel((base + std::to_string(iCall++) + ".lp").c_str());

  cplex.setOut(env.getNullStream());
  if ( !cplex.solve() ) {
    env.error() << "Failed to optimize LP" << endl;
    exit(-1);
    return false;
  }
  
  bool verbose = false;

  IloNumArray vals(env);
  if (verbose) env.out() << "Solution status = " << cplex.getStatus() << endl;
  if (verbose) env.out() << "Solution value  = " << cplex.getObjValue() << endl;
  // cplex.getValues(vals, cols);
  //env.out() << "Values        = " << vals << endl;

  for (int i = 0; i < nverts; i++) {
    double x[3] = { 0, 0, 0 };
    for (int xyz = 0; xyz < 3; xyz++)
      for (int pm = 0; pm < 2; pm++)
        // x[xyz] += vals[6*i + 2*xyz + pm] * (1-2*pm);
	x[xyz] += cplex.getValue(cols[6*i + 2*xyz + pm]) * (1-2*pm);
    values[verts[i]] = Point(x);

    if (verbose)
      cerr << verts[i] << " "
	   << values[verts[i]].x[0] << " "
	   << values[verts[i]].x[1] << " "
	   << values[verts[i]].x[2] << endl;
  }

  env.end();
  
  return true;
}

double expanderT;

double saveV;
vector<bool> skipPair;

bool Expander2::expandV (double e, bool velocityObjective, double velocityBound) {
  if (true) {
    double uscale = uscaleBase;
    double violation;
    while ((violation = expandV2(e, velocityObjective, velocityBound, uscale)) > 1e-6) {
      if (violation == 1234567890) {
	velocityBound /= 2;
	// cerr << "Feature intersection.  Restarting with velocityBound " << velocityBound << endl;
	uscale *= 10;
	cerr << "Feature intersection. Restarting with uscale " << uscale << endl;
      }
      else {
	uscale *= 10;
	cerr << "row violation " << violation << " restarting with uscale " << uscale << endl;
      }
    }
    return true;
  }

  if (false) {
    skipPair.resize(pairs.size());
    for (int i = 0; i < skipPair.size(); i++)
      skipPair[i] = false;
    
    for (int i = 0; i < skipPair.size(); i++) {
      skipPair[i] = true;
      expandV2(e, velocityObjective, velocityBound, 1);
      if (saveV > 0)
	skipPair[i] = false;
    }
    for (int i = 0; i < skipPair.size(); i++)
      if (!skipPair[i])
	cerr << "pair " << i << " " << pairs[i].d << endl;
  }
}

double minSep1;

double Expander2::expandV2 (double e, bool velocityObjective, double velocityBound, double uscale) {
  bool verbose = false;

  minSep1 = 1e9;

  values.clear();
  assert(0 < e && e <= 1);
  assert(0 < velocityBound && velocityBound <= 1);

  static int iCall = 100;

  vector<Vertex*> verts;
  map<Vertex*, int> index;

  int npairs = pairs.size();

  if (verbose)
    cerr << "vobj " << velocityObjective << " vbnd " << velocityBound << endl;

  int isep = -1;
  double minsep = 100;
  for (int ipair = 0; ipair < npairs; ipair++) {
    if (pairs[ipair].d < minsep) {
      isep = ipair;
      minsep = pairs[ipair].d;
    }
  }
  if (verbose) {
    cerr << iCall << " minimum separation " << isep << " " << minsep << endl;
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[isep].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
	cerr << constraints[ivert].i << " ";
      }
      cerr << endl;
    }
  }

  /*
  cerr << "476 pairs:" << endl;
  for (int ipair = 0; ipair < npairs; ipair++)
    if (pairs[ipair].constraints[0][0].i == 476) {
      for (int ifeat = 0; ifeat < 2; ifeat++) {
	vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
	for (int ivert = 0; ivert < constraints.size(); ivert++) {
	  cerr << constraints[ivert].i << " ";
	}
	cerr << endl;
      }
    }
  */

  double totalDisp = 0;

  for (int ipair = 0; ipair < npairs; ipair++)
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Vertex* vert = constraints[ivert].i;
        if (index.find(vert) == index.end()) {
          index[vert] = verts.size();
          verts.push_back(vert);

	  Point disp = disps[vert];
	  totalDisp += fabs(disp.x[0]) + fabs(disp.x[1]) + fabs(disp.x[2]);
        }
      }
    }

  int nverts = verts.size();

  if (verbose)
    cerr << "npairs " << npairs << " nverts " << nverts << " displacement " << totalDisp << endl;

#ifdef BLEEN
  static double saveDisp = 1e6;
  static double saveDisp2;
  static double minimizeT = 1;
  // bool projection = minsep < 0.99999;

  bool projection = (saveDisp == 1e6 && minsep < 0.99999 || iCall % 2 == 0);
  if (!projection) {
    double ratio = (saveDisp - totalDisp) / (saveDisp - saveDisp2);
    /*
    static int rcount;
    if (ratio < 0.25) {
      minimizeT /= 2;
      rcount = 0;
    }
    else if (ratio > 0.5) {
      rcount++;
      if (rcount >= 4 && minimizeT < 1) {
	minimizeT *= 2;
	rcount = 0;
      }
    }
    */

    if (ratio > 0.5) {
      minimizeT *= 1.5;
      if (minimizeT > 1)
	minimizeT = 1;
    }
    else
      minimizeT /= 1.5;

    /*
    ifstream in;
    in.open("minimizet");
    in >> minimizeT;
    in.close();
    */

    saveDisp = totalDisp;
    // expanderT = minimizeT;
    //cerr << "ratio " << ratio << " minimize t " << minimizeT << endl;
  }
  else {
    expanderT = 1;
    saveDisp2 = totalDisp;
  }
#endif

  IloEnv env;
  IloNumVarArray cols(env);
  IloNumArray coefs(env);

  // a' or b' is (x,y,z)
  // -|a| <= a + a' < |a|
  string xyzll[] = { "x", "lxl", "y", "lyl", "z", "lzl" };
  for (int i = 0; i < nverts; i++) {
    string base = "v";
    base += std::to_string(i);
    for (int j = 0; j < 6; j++) {
      if (j % 2 == 0)
        // cols.add(IloNumVar(env, -1, 1));
	cols.add(IloNumVar(env, -velocityBound, velocityBound));
      else
        cols.add(IloNumVar(env));
      cols[6*i+j].setName((base + xyzll[j]).c_str());
      if (j % 2 == 0)
        coefs.add(IloNum(0));
      else
        coefs.add(IloNum(1));
    }
  }

  // p'*u, q'*u, and u' = v' v + w' w
  string pvw[] = { "p", "q", "v", "w" };
  for (int i = 0; i < npairs; i++) {
    Pair &pair = pairs[i];
    string base = "p";
    base += std::to_string(i);
    for (int j = 0; j < 4; j++) {
      if (j < 2)
	cols.add(IloNumVar(env, -sqrt(3.0), sqrt(3.0)));
      else if (j == 2)
	cols.add(IloNumVar(env, pair.lMin, pair.lMax));
      else if (j == 3)
	cols.add(IloNumVar(env, pair.mMin, pair.mMax));
      cols[6*nverts + 4*i + j].setName((base + pvw[j]).c_str());
      coefs.add(IloNum(0));
    }
  }

  int VCol = 6*nverts + 4*npairs;
  cols.add(IloNumVar(env, 0, 1));
  cols[VCol].setName("V");
  // coefs.add(IloNum(-6 * nverts));
  coefs.add(IloNum(-6000 * nverts));

  IloObjective obj = IloMinimize(env);

  obj.setLinearCoefs(cols, coefs);
    
  IloRangeArray rows(env);
  vector<Ineq> ineqs;

  int nrows = 0;
  for (int ipair = 0; ipair < npairs; ipair++) {
    // if (skipPair[ipair])
    // continue;

    string base = "p";
    base += std::to_string(ipair);
    Pair &pair = pairs[ipair];
    int pCol = 6 * nverts + 4 * ipair;
    
    if (false) {
      // qu' - pu' >=  (1 - (q - p) * u / d) * V
      // qu' - pu' - (e - (q - p) * u / d) * V >= 0
      rows.add(IloRange(env, 0, IloInfinity));
      ineqs.push_back(Ineq(true, 0));
      IloRange row = rows[nrows];
      nrows++;
      row.setName((base + "pq").c_str());
      row.setLinearCoef(cols[pCol], -1);
      ineqs.back().add(pCol, -1);
      row.setLinearCoef(cols[pCol+1], 1);
      ineqs.back().add(pCol+1, 1);
      row.setLinearCoef(cols[VCol], -(e - pair.d));
      ineqs.back().add(VCol, -(e - pair.d));
    }
    else if (false) {
      // qu' - pu' >=  e * V - (q - p) * u / d
      // qu' - pu' - e * V >= - (q - p) * u / d
      rows.add(IloRange(env, -pair.d, IloInfinity));
      ineqs.push_back(Ineq(true, -pair.d));
      IloRange row = rows[nrows];
      nrows++;
      row.setName((base + "pq").c_str());
      row.setLinearCoef(cols[pCol], -1);
      ineqs.back().add(pCol, -1);
      row.setLinearCoef(cols[pCol+1], 1);
      ineqs.back().add(pCol+1, 1);
      row.setLinearCoef(cols[VCol], -e);
      ineqs.back().add(VCol, -e);
    }

    for (int ifeat = 0; ifeat < 0; ifeat++) {
      vector<Constraint> &constraints = pair.constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Constraint con = constraints[ivert];

	// if (con.r != 0)
	// continue;

        int vCol = 6 * index[con.i];
        if (ifeat == 0) {
          rows.add(IloRange(env, -IloInfinity, -con.r));
          ineqs.push_back(Ineq(false, -con.r));
          rows[nrows].setName((base + "a" + std::to_string(index[con.i])).c_str());
        }
        else {
          rows.add(IloRange(env, -con.r, IloInfinity));
          ineqs.push_back(Ineq(true, -con.r));
          rows[nrows].setName((base + "b" + std::to_string(index[con.i])).c_str());
        }
        
        IloRange row = rows[nrows];

        // a' * u - pu' + x' (a - p) * v / d + y' (a - p) * w / d <= (-(a - p) * u / d) * V
        // a' * u - pu' + x' (a - p) * v / d + y' (a - p) * w / d + V (a - p) * u / d
        // (a' - p') * u + x (a - p) * v + y (a - p) * w
        for (int xyz = 0; xyz < 3; xyz++) {
          row.setLinearCoef(cols[vCol + 2 * xyz], pair.u.x[xyz]);
          ineqs.back().add(vCol + 2 * xyz, pair.u.x[xyz]);
        }


	row.setLinearCoef(cols[pCol + ifeat], -1);
        ineqs.back().add(pCol + ifeat, -1);

	if (true) {
	  // double uscale = 1;
	  row.setLinearCoef(cols[pCol + 2],  con.v / uscale);
	  ineqs.back().add(pCol + 2,  con.v / uscale);
	  if (con.r == 0) {
	    row.setLinearCoef(cols[pCol + 3],  con.w / uscale);
	    ineqs.back().add(pCol + 3,  con.w / uscale);
	  }
	}

        // row.setLinearCoef(cols[VCol], con.r);
        // ineqs.back().add(VCol, con.r);

        nrows++;
      }
    }

    // number of vertices on separating planes (r=0)
    int numOnSep = 0;
    for (int ivertA = 0; ivertA < pair.constraints[0].size(); ivertA++) {
      Constraint conA = pair.constraints[0][ivertA];
      if (conA.r == 0)
        numOnSep++;
    }
    for (int ivertB = 0; ivertB < pair.constraints[1].size(); ivertB++) {
      Constraint conB = pair.constraints[1][ivertB];
      if (conB.r == 0)
        numOnSep++;
    }

    for (int ivertA = 0; ivertA < pair.constraints[0].size(); ivertA++) {
      Constraint conA = pair.constraints[0][ivertA];
      for (int ivertB = 0; ivertB < pair.constraints[1].size(); ivertB++) {
        Constraint conB = pair.constraints[1][ivertB];

        int vColA = 6 * index[conA.i];
        int vColB = 6 * index[conB.i];

        rows.add(IloRange(env, -pair.d + conA.r - conB.r, IloInfinity));
        ineqs.push_back(Ineq(true, -pair.d + conA.r - conB.r));
        rows[nrows].setName((base + 
                             "a" + std::to_string(index[conA.i]) +
                             "b" + std::to_string(index[conB.i])).c_str());
        IloRange row = rows[nrows];

        row.setLinearCoef(cols[VCol], -e);
        ineqs.back().add(VCol, -e);

        for (int xyz = 0; xyz < 3; xyz++) {
          row.setLinearCoef(cols[vColA + 2 * xyz], -pair.u.x[xyz]);
          ineqs.back().add(vColA + 2 * xyz, -pair.u.x[xyz]);
        }

        for (int xyz = 0; xyz < 3; xyz++) {
          row.setLinearCoef(cols[vColB + 2 * xyz], pair.u.x[xyz]);
          ineqs.back().add(vColB + 2 * xyz, pair.u.x[xyz]);
        }

        if (!pair.limitVirtualPairU || numOnSep >= 3) {
          row.setLinearCoef(cols[pCol + 2],  (conB.v - conA.v) / uscale);
          ineqs.back().add(pCol + 2,  (conB.v - conA.v) / uscale);
        }
        if (!pair.limitVirtualPairU || numOnSep == 4) {
          row.setLinearCoef(cols[pCol + 3],  (conB.w - conA.w) / uscale);
          ineqs.back().add(pCol + 3,  (conB.w - conA.w) / uscale);
        }

	nrows++;
      }
    }
  }

  for (int ivert = 0; ivert < nverts; ivert++) {
    string base = "v";
    base += std::to_string(ivert);
    Point disp = disps[verts[ivert]];
    int vCol = 6 * ivert;

    string xyz[] = { "x", "y", "z" };
    for (int i = 0; i < 3; i++) {
      int xCol = vCol + 2 * i;
      int lxlCol = xCol + 1;
      // double d = iCall % 2 == 0 ? disp.x[i] : 0;
      double d = velocityObjective ? 0 : disp.x[i];

      // d + x < |x|
      // -x + |x| > d
      rows.add(IloRange(env, d, IloInfinity));
      ineqs.push_back(Ineq(true, d));
      IloRange row = rows[nrows];
      nrows++;
      row.setName((base + xyz[i] + "u").c_str());
      row.setLinearCoef(cols[xCol], -1.0);
      ineqs.back().add(xCol, -1.0);
      row.setLinearCoef(cols[lxlCol], 1.0);
      ineqs.back().add(lxlCol, 1.0);

      // d + x > -|x|
      // x + |x| > -d
      rows.add(IloRange(env, -d, IloInfinity));
      ineqs.push_back(Ineq(true, -d));
      row = rows[nrows];
      nrows++;
      row.setName((base + xyz[i] + "l").c_str());
      row.setLinearCoef(cols[xCol], 1.0);
      ineqs.back().add(xCol, 1.0);
      row.setLinearCoef(cols[lxlCol], 1.0);
      ineqs.back().add(lxlCol, 1.0);
    }
  }

  IloModel model(env);

  model.add(obj);
  model.add(rows);

  IloCplex cplex(model);

  string base = "expand";
  // cplex.exportModel((base + std::to_string(iCall++) + ".lp").c_str());

  cplex.setOut(env.getNullStream());
  //cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
  if ( !cplex.solve() ) {
    env.error() << "Failed to optimize LP" << endl;
    // exit(-1);
    return 1;
  }
  
  bool verbose2 = false;

  IloNumArray vals(env);
  if (verbose) env.out() << "Solution status = " << cplex.getStatus() << endl;
  if (verbose) env.out() << "Solution value  = " << cplex.getObjValue() << endl;
  // cplex.getValues(vals, cols);
  //env.out() << "Values        = " << vals << endl;

  double maxval = 0;
  int imaxval = -1;
  for (int i = 0; i < ineqs.size(); i++) {
    double v = ineqs[i].value(cplex, cols);
    if (v > maxval) {
      imaxval = i;
      maxval = v;
    }
  }

  if (verbose && imaxval >= 0 && (verbose || maxval > 1))
    cerr << (iCall-1) << " row violation " << imaxval << " " << rows[imaxval].getName() << " " << maxval << endl;

  if (verbose)
    cerr << "V " << cplex.getValue(cols[VCol]) << endl;
  
  double maxT = 1;
  for (int i = 0; i < nverts; i++) {
    double x[3] = { 0, 0, 0 };
    for (int xyz = 0; xyz < 3; xyz++)
      for (int pm = 0; pm < 1; pm++)
        // x[xyz] += vals[6*i + 2*xyz + pm] * (1-2*pm);
	x[xyz] += cplex.getValue(cols[6*i + 2*xyz + pm]) * (1-2*pm) * maxT;
    values[verts[i]] = Point(x);

    if (verbose2 && fabs(x[0]) + fabs(x[1]) + fabs(x[2]) > 1e-8)
      cerr << verts[i] << " "
	   << values[verts[i]].x[0] << " "
	   << values[verts[i]].x[1] << " "
	   << values[verts[i]].x[2] << endl;
  }

  if (maxval < 1e-6) {
    for (int i = 0; i < pairs.size(); i++) {
      while (!checkPair(cplex, cols, index, i, maxT, uscale))
        maxT *= 0.9;
    }
    // maxT = checkPair2(cplex, cols, index, i, maxT, uscale);
    if (verbose)
      cerr << "maxT " << maxT << endl;
  }
  if (verbose)
    cerr << endl;

  if (maxT < 1.0) {
    for (int i = 0; i < nverts; i++) {
      double x[3] = { 0, 0, 0 };
      values[verts[i]] = Point(x);
    }
    return 1234567890;
  }

  // maxT = 1;

  //  saveV = cplex.getValue(cols[VCol] * maxT);

  env.end();
  
  return maxval;
}

bool Expander2::checkPair (IloCplex &cplex, IloNumVarArray& cols,
                           map<Vertex*, int> &index,
                           int ipair, double t, double s) {
  bool verbose = false;
  Pair &pair = pairs[ipair];
  int pCol = 6 * index.size() + 4 * ipair;

  double p_ = 0; //cplex.getValue(cols[pCol]);
  double q_ = 0; //cplex.getValue(cols[pCol+1]);

  double sep = pair.d + t * (q_ - p_);

  double sep1 = pair.d + (q_ - p_);
  if (minSep1 > sep1) {
    if (verbose)
      cerr << "minSep1 " << sep1 << endl;
    minSep1 = sep1;
  }

  double x_ = 0, y_ = 0;
  for (int ab = 0; ab <= 1; ab++) {
    vector<Constraint> &constraints = pair.constraints[ab];
    for (int i = 0; i < constraints.size(); i++) {
      Constraint &c = constraints[i];
      try {
	x_ = c.v != 0 ? cplex.getValue(cols[pCol+2]) : x_;
	y_ = c.r == 0 && c.w != 0 ? cplex.getValue(cols[pCol+3]) : y_;
      } catch (IloAlgorithm::NotExtractedException nee) {
      }
    }
  }

#ifdef BLEEN
  double disMax[2] = { 0, 0 };
  for (int ab = 0; ab <= 1; ab++) {
    vector<Constraint> &constraints = pair.constraints[ab];
    for (int i = 0; i < constraints.size(); i++) {
      Constraint &c = constraints[i];
      int vCol = 6 * index[c.i];
      Point a_(cplex.getValue(cols[vCol+0]),
               cplex.getValue(cols[vCol+2]),
               cplex.getValue(cols[vCol+4]));
      double pq_ = ab == 0 ? p_ : q_;
      double dis = c.r +
        t * (x_ * c.v / s + y_ * c.w / s + a_.dot(pair.u) - pq_) +
        t * t * (x_ * a_.dot(pair.v) / s + y_ * a_.dot(pair.w) / s);
   
      double dis1 = c.r + x_ * c.v / s + y_ * c.w / s + a_.dot(pair.u) - pq_;
      if (ab == 0) {
        if (i == 0 || disMax[0] < dis)
          disMax[0] = dis;

	static double maxDis1 = -1e9;
	if (maxDis1 < dis1) {
	  if (verbose)
	    cerr << "maxDis1 " << dis1 << endl;
	  maxDis1 = dis1;
	}
      }
      else {
        if (i == 0 || disMax[1] > dis)
          disMax[1] = dis;

	static double minDis1 = 1e9;
	if (minDis1 > dis1) {
	  if (verbose)
	    cerr << "minDis1 " << dis1 << endl;
	  minDis1 = dis1;
	}
      }
    }
  }
#endif

  double mindist = 1e9;
  for (int iA = 0; iA < pair.constraints[0].size(); iA++) {
      Constraint &cA = pair.constraints[0][iA];
      int vColA = 6 * index[cA.i];
      Point a_(cplex.getValue(cols[vColA+0]),
               cplex.getValue(cols[vColA+2]),
               cplex.getValue(cols[vColA+4]));
      for (int iB = 0; iB < pair.constraints[1].size(); iB++) {
        Constraint &cB = pair.constraints[1][iB];
        int vColB = 6 * index[cB.i];
        Point b_(cplex.getValue(cols[vColB+0]),
                 cplex.getValue(cols[vColB+2]),
                 cplex.getValue(cols[vColB+4]));
      double dist = pair.d + cB.r - cA.r +
        t * (x_ * (cB.v - cA.v) / s + y_ * (cB.w - cA.w) / s + b_.dot(pair.u) - a_.dot(pair.u)) +
        t * t * (x_ * (b_.dot(pair.v) - a_.dot(pair.v)) / s + 
                 y_ * (b_.dot(pair.w) - a_.dot(pair.w)) / s);
   
      if (dist < mindist)
        mindist = dist;
    }
  }

  double sep2 = mindist; // sep - disMax[0] + disMax[1]; // 
  static double minSep = 1e9;
  if (minSep > sep2) {
    if (verbose)
      cerr << "minSep " << sep2 << " d " << pair.d << " sep1 " << sep1 << endl;
    minSep = sep2;
  }
  if (sep2 < 0)
    minSep = 1e9;
  bool intersects (Expander2 *, Expander2::Pair *);
  return sep2 > 0 || !intersects(this, &pair);
}

#include "pv.h"
using namespace acp;

double Expander2::checkPair2 (IloCplex &cplex, IloNumVarArray& cols,
			      map<Vertex*, int> &index,
			      int ipair, double t, double s) {
  bool verbose = false;
  Pair &pair = pairs[ipair];

  vector<Point> ps[2];
  vector<Point> vs[2];
  int abs[4], is[4], n = 0;
  for (int ab = 0; ab <= 1; ab++) {
    vector<Constraint> &constraints = pair.constraints[ab];
    for (int i = 0; i < constraints.size(); i++) {
      assert(n < 4);
      abs[n] = ab;
      is[n] = i;
      n++;
      
      Constraint &c = constraints[i];
      if (fabs(c.r) > 0.1)
	return t;
      
      double x = c.v;
      double y = c.w;
      double z = ab ? c.r + pair.d : c.r;
      ps[ab].push_back(Point(x, y, z));

      int vCol = 6 * index[c.i];
      Point a_(Point(cplex.getValue(cols[vCol+0]),
		     cplex.getValue(cols[vCol+2]),
		     cplex.getValue(cols[vCol+4])));
      vs[ab].push_back(Point(pair.v.dot(a_), pair.w.dot(a_), pair.u.dot(a_)));
    }
  }
  assert(n == 0);

  Point p01 = ps[abs[1]][is[1]] - ps[abs[0]][is[0]];
  Point p02 = ps[abs[2]][is[2]] - ps[abs[0]][is[0]];
  Point p03 = ps[abs[3]][is[3]] - ps[abs[0]][is[0]];
  Point v01 = vs[abs[1]][is[1]] - vs[abs[0]][is[0]];
  Point v02 = vs[abs[2]][is[2]] - vs[abs[0]][is[0]];
  Point v03 = vs[abs[3]][is[3]] - vs[abs[0]][is[0]];
  
  double cs[4];
  cs[0] = (p01.cross(p02).dot(p03));
  cs[1] = (v01.cross(p02).dot(p03) +
	   p01.cross(v02).dot(p03) +
	   p01.cross(p02).dot(v03));
  cs[2] = (p01.cross(v02).dot(v03) +
	   v01.cross(p02).dot(v03) +
	   v01.cross(v02).dot(p03));
  cs[3] = (v01.cross(v02).dot(v03));

  double C = cs[1];
  double B = cs[2] * 2;
  double A = cs[3] * 3;
  if (A < 0) {
    C = -C;
    B = -B;
    A = -A;
  }
  
  double D = B * B - A * C * 4;

  if (D < 0)
    return t;

  double r[2];
  if (B > 0) {
    double BD = -B - sqrt(D);
    r[0] = BD / (A * 2);
    r[1] = (C * 2) / BD;
  }
  else {
    double BD = -B + sqrt(D);
    r[0] = (C * 2) / BD;
    r[1] = BD / (A * 2);
  }
  
  double saveT = t;
  
  if (r[0] > 0 && r[0] < t)
    t = r[0];
  else if (r[1] > 0 && r[1] < t)
    t = r[1];
  
  if (verbose && t < saveT)
    cerr << "t " << t << endl;
    
  return t;
}
