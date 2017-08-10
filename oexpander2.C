#include "expander2.h"
#include <ilcplex/ilocplex.h>
#include <string>

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
      // cout << "extracting " << cols[i] << " " << ccols[cols[i]].getName() << endl;
      v += cplex.getValue(ccols[cols[i]]) * coefs[i];
    }
    // cout << endl;
    return gt ? -v : v;
  }
};

bool Expander2::expand () {
  vector<int> verts;
  map<int, int> index;

  int npairs = pairs.size();

  for (int ipair = 0; ipair < npairs; ipair++)
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        int vert = constraints[ivert].i;
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
    base += std::to_string(verts[i]);
    for (int j = 0; j < 6; j++) {
      cols.add(IloNumVar(env));
      // cout << "cols.size() " << cols.length << " = " << 6*i+j << endl;
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
          rows[nrows].setName((base + "a" + std::to_string(con.i)).c_str());
        }
        else {
          rows.add(IloRange(env, 1 - con.r, IloInfinity));
          rows[nrows].setName((base + "b" + std::to_string(con.i)).c_str());
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

  cplex.exportModel("expand.lp");

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
      cout << verts[i] << " "
	   << values[verts[i]].x[0] << " "
	   << values[verts[i]].x[1] << " "
	   << values[verts[i]].x[2] << endl;
  }

  env.end();
  
  return true;
}

bool Expander2::expand (double e) {
  vector<int> verts;
  map<int, int> index;

  int npairs = pairs.size();

  for (int ipair = 0; ipair < npairs; ipair++)
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        int vert = constraints[ivert].i;
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
    base += std::to_string(verts[i]);
    for (int j = 0; j < 6; j++) {
      cols.add(IloNumVar(env, 0, 1));
      // cout << "cols.size() " << cols.length << " = " << 6*i+j << endl;
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
    cout << "setting objective coefficients of vertices" << endl;
    for (int i = 0; i < nverts; i++)
      for (int j = 0; j < 6; j++)
	obj.setLinearCoef(cols[6*i+j], 1.0);
    cout << "done" << endl;
    
    cout << "setting objective coefficients of pairs" << endl;
    bool uObj = true;
    if (uObj) {
      for (int i = 0; i < npairs; i++)
	for (int j = 0; j < 6; j++)
	  obj.setLinearCoef(cols[6*nverts+6*i+j], 1.0);
    }
    cout << "done" << endl;
    
    obj.setLinearCoef(cols[dCol], -6 * nverts);
  }
  else {
    //cout << "setting objective" << endl;
    obj.setLinearCoefs(cols, coefs);
    //cout << "done" << endl;
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
          rows[nrows].setName((base + "a" + std::to_string(con.i)).c_str());
        }
        else {
          rows.add(IloRange(env, - con.r, IloInfinity));
          rows[nrows].setName((base + "b" + std::to_string(con.i)).c_str());
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
  cplex.exportModel((base + std::to_string(iCall++) + ".lp").c_str());

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
      cout << verts[i] << " "
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
      uscale *= 10;
      cout << "row violation " << violation << " restarting with uscale " << uscale << endl;
    }
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
	cout << "pair " << i << " " << pairs[i].d << endl;
  }
}


double Expander2::expandV2 (double e, bool velocityObjective, double velocityBound, double uscale) {
  values.clear();
  bool verbose = false;
  assert(0 < e && e <= 1);
  assert(0 < velocityBound && velocityBound <= 1);

  static int iCall = 100;

  vector<int> verts;
  map<int, int> index;

  int npairs = pairs.size();

  if (verbose)
    cout << "vobj " << velocityObjective << " vbnd " << velocityBound << endl;

  int isep = -1;
  double minsep = 100;
  for (int ipair = 0; ipair < npairs; ipair++) {
    if (pairs[ipair].d < minsep) {
      isep = ipair;
      minsep = pairs[ipair].d;
    }
  }
  if (verbose) {
    cout << iCall << " minimum separation " << isep << " " << minsep << endl;
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[isep].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
	cout << constraints[ivert].i << " ";
      }
      cout << endl;
    }
  }

  /*
  cout << "476 pairs:" << endl;
  for (int ipair = 0; ipair < npairs; ipair++)
    if (pairs[ipair].constraints[0][0].i == 476) {
      for (int ifeat = 0; ifeat < 2; ifeat++) {
	vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
	for (int ivert = 0; ivert < constraints.size(); ivert++) {
	  cout << constraints[ivert].i << " ";
	}
	cout << endl;
      }
    }
  */

  double totalDisp = 0;

  for (int ipair = 0; ipair < npairs; ipair++)
    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pairs[ipair].constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        int vert = constraints[ivert].i;
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
    cout << "npairs " << npairs << " nverts " << nverts << " displacement " << totalDisp << endl;

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
    //cout << "ratio " << ratio << " minimize t " << minimizeT << endl;
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
    base += std::to_string(verts[i]);
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
    string base = "p";
    base += std::to_string(i);
    for (int j = 0; j < 4; j++) {
      if (j < 2)
	cols.add(IloNumVar(env, -sqrt(3.0), sqrt(3.0)));
      else
	cols.add(IloNumVar(env, -1, 1));
      cols[6*nverts + 4*i + j].setName((base + pvw[j]).c_str());
      coefs.add(IloNum(0));
    }
  }

  int VCol = 6*nverts + 4*npairs;
  cols.add(IloNumVar(env, 0, 1));
  cols[VCol].setName("V");
  // coefs.add(IloNum(-6 * nverts));
  coefs.add(IloNum(-60 * nverts));

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

    for (int ifeat = 0; ifeat < 2; ifeat++) {
      vector<Constraint> &constraints = pair.constraints[ifeat];
      for (int ivert = 0; ivert < constraints.size(); ivert++) {
        Constraint con = constraints[ivert];

	// if (con.r != 0)
	// continue;

        int vCol = 6 * index[con.i];
        if (ifeat == 0) {
          rows.add(IloRange(env, -IloInfinity, -con.r));
          ineqs.push_back(Ineq(false, -con.r));
          rows[nrows].setName((base + "a" + std::to_string(con.i)).c_str());
        }
        else {
          rows.add(IloRange(env, -con.r, IloInfinity));
          ineqs.push_back(Ineq(true, -con.r));
          rows[nrows].setName((base + "b" + std::to_string(con.i)).c_str());
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
	  row.setLinearCoef(cols[pCol + 3],  con.w / uscale);
	  ineqs.back().add(pCol + 3,  con.w / uscale);
	}

        // row.setLinearCoef(cols[VCol], con.r);
        // ineqs.back().add(VCol, con.r);

        nrows++;
      }
    }
  }

  for (int ivert = 0; ivert < nverts; ivert++) {
    string base = "v";
    base += std::to_string(verts[ivert]);
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
  cplex.exportModel((base + std::to_string(iCall++) + ".lp").c_str());

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

  if (imaxval >= 0 && (verbose || maxval > 1))
    cout << (iCall-1) << " row violation " << imaxval << " " << rows[imaxval].getName() << " " << maxval << endl;

  if (verbose)
    cout << "V " << cplex.getValue(cols[VCol]) << endl << endl;
  
  saveV = cplex.getValue(cols[VCol]);

  for (int i = 0; i < nverts; i++) {
    double x[3] = { 0, 0, 0 };
    for (int xyz = 0; xyz < 3; xyz++)
      for (int pm = 0; pm < 1; pm++)
        // x[xyz] += vals[6*i + 2*xyz + pm] * (1-2*pm);
	x[xyz] += cplex.getValue(cols[6*i + 2*xyz + pm]) * (1-2*pm);
    values[verts[i]] = Point(x);

    if (verbose2 && fabs(x[0]) + fabs(x[1]) + fabs(x[2]) > 1e-8)
      cout << verts[i] << " "
	   << values[verts[i]].x[0] << " "
	   << values[verts[i]].x[1] << " "
	   << values[verts[i]].x[2] << endl;
  }

  env.end();
  
  return maxval;
}
