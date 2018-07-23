#include "sumfaces.h"
#include "io.h"
#include <time.h>
#include <string>

const bool BLOCKSPACES_FROM_FILE = true;
const bool GENERATE_INPUT_FILE = false;
const bool GENERATE_ANGLE_RANGES_FILE = false;
const bool COMPARE_LIFESPANS = true;

double tanHalfAngle (int n) {
  return tan((1.0 + 1.0e-8) * M_PI / n);
}

double sampleTanHalfAngle(int num_rotations, int samples_per_rotation) {
  double sample_theta = (2 * M_PI / num_rotations - 1e-8) / samples_per_rotation;
  return tan(sample_theta / 2);
}

void generateInputFile(const char * directory) {
  std::string dir(directory);

  Polyhedron * robot = loadPoly((dir + "/robot").c_str());
  Polyhedron * obstacle = loadPoly((dir + "/obstacle").c_str());

  ofstream out;
  out.open((dir + "/pairs.txt").c_str());
  for (int i=0; i<robot->vertices.size(); i++)
    for (int j=0; j<obstacle->faces.size(); j++) {
      out << "3 ";
      Vertices vs;
      obstacle->faces[j]->boundaryVertices(vs);
      assert(vs.size() == 3);
      for (int k=0; k<3; k++)
        out << (int)(std::find(obstacle->vertices.begin(), obstacle->vertices.end(), vs[k]) - obstacle->vertices.begin()) << " ";
      out<< "1 " << i << endl;
    }

  for (int i=0; i<robot->faces.size(); i++)
    for (int j=0; j<obstacle->vertices.size(); j++) {
      out << "1 "<< j << " 3 ";
      Vertices vs;
      robot->faces[i]->boundaryVertices(vs);
      assert(vs.size() == 3);
      for (int k=0; k<3; k++)
        out << (int)(std::find(robot->vertices.begin(), robot->vertices.end(), vs[k]) - robot->vertices.begin()) << " ";
      out << endl;
    }

  for (int i=0; i<robot->edges.size(); i++)
    for (int j=0; j<obstacle->edges.size(); j++) {
      out << "2 ";
      out << (int)(find(obstacle->vertices.begin(), obstacle->vertices.end(), obstacle->edges[j]->getH()) - obstacle->vertices.begin()) << " ";
      out << (int)(find(obstacle->vertices.begin(), obstacle->vertices.end(), obstacle->edges[j]->getT()) - obstacle->vertices.begin()) << " ";
      out << "2 ";
      out << (int)(find(robot->vertices.begin(), robot->vertices.end(), robot->edges[i]->getH()) - robot->vertices.begin()) << " ";
      out << (int)(find(robot->vertices.begin(), robot->vertices.end(), robot->edges[i]->getT()) - robot->vertices.begin()) << endl;
    }

  out.close();
} 

double angle(double x, double y) {
  Parameter::disable();
  double theta;

  if (x > 1.0)
    x = 1.0;
  if (x < -1.0)
    x = -1.0;
  if (y > 1.0)
    y = 1.0;
  if (y < -1.0)
    y = -1.0;

  if (x >= 0 && y >= 0)
    theta = acos(x);
  else if (x <= 0 && y >= 0)
    theta = M_PI/2 + acos(y);
  else if (x <= 0 && y <=0)
    theta = M_PI + acos(-x);
  else
    theta = 2 * M_PI - acos(x);

  Parameter::enable();

  return theta;
}

bool * lifespansToAngleRanges(std::vector<std::pair<double, double> > & lifespans, int num_rotations) {
  double theta = 2*(1.0 + 1.0e-8) * M_PI / num_rotations;
  bool * candidate =  new bool[num_rotations];
  for (int i=0; i<num_rotations; i++)
    candidate[i] = false;

  for (int i=0; i<lifespans.size(); i++) {
    int start_range = floor(lifespans[i].first / theta);
    int end_range = floor(lifespans[i].second / theta);

    if (lifespans[i].first <= lifespans[i].second) {
      for (int j=start_range; j<=end_range; j++)
        candidate[j] = true;
    } else {
      for (int j=start_range; j<num_rotations; j++)
        candidate[j] = true;
      for (int j=0; j<end_range; j++)
        candidate[j] = true;
    }
  }

  return candidate;
}

std::string getKey(int n[6]) {
  int tmp;
  if (n[0] == 1) {
    //sort 345
    if (n[3] > n[4]) {
      tmp = n[3];
      n[3] = n[4];
      n[4] = tmp;
    }
    if (n[3] > n[5]) {
      tmp = n[3];
      n[3] = n[5];
      n[5] = tmp;
    }
    if (n[4] > n[5]) {
      tmp = n[4];
      n[4] = n[5];
      n[5] = tmp;
    }
  } else if (n[0] == 2) {
    if (n[1] > n[2]) {
      tmp = n[1];
      n[1] = n[2];
      n[2] = tmp;
    }
    if (n[4] > n[5]) {
      tmp = n[4];
      n[4] = n[5];
      n[5] = tmp;
    }
  } else if (n[0] == 3) {
    //sort 123
    if (n[1] > n[2]) {
      tmp = n[1];
      n[1] = n[2];
      n[2] = tmp;
    }
    if (n[1] > n[3]) {
      tmp = n[1];
      n[1] = n[3];
      n[3] = tmp;
    }
    if (n[2] > n[3]) {
      tmp = n[2];
      n[2] = n[3];
      n[3] = tmp;
    } 
  } else {
    cout << "invalid pair" << endl;
    return NULL;
  }

  return std::to_string(n[0])+" "+std::to_string(n[1])+" "+std::to_string(n[2])+" "+std::to_string(n[3])+" "+std::to_string(n[4])+" "+std::to_string(n[5]);
}

void createMap(const char * pairfile, const char * outfile, std::map<std::string, bool*> & pairmap) {
  ifstream infile1(pairfile);
  ifstream infile2(outfile);

  std::string s;
  int n[6];
  while (infile1 >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5]) {
    infile2 >> s;
    assert(s.length() == 40);

    bool * b = new bool[s.length()];
    for (int i=0; i<s.length(); i++) {
      if (s[i] == '0')
        b[i] = false;
      else if (s[i] == '1')
        b[i] = true;
      else 
        assert(false);
    }

    std::string key = getKey(n);
    pairmap[key] = b;
    cout << "added key '"<<key.c_str()<<"'"<<endl;
  }
}

void testLifeSpans(std::string directory, std::string outfile, int num_rotations) {
  std::map<std::string, bool *> pairmap;
  createMap((directory+"/pairs.txt").c_str(), (directory+"/"+outfile).c_str(), pairmap);

  ifstream infile ((directory+"/faces.txt").c_str());
  int n[6];
  char s[50];
  int lifespan_candidates = 0;
  int eliminated_candidates = 0;
    

  while (infile >> s >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5]) {
    assert(strcmp(s,"face") == 0);
    cout << s << ": " << n[0] << " " << n[1] << " " << n[2] << " " << n[3] << " " << n[4] << " " << n[5] << endl;

    std::string key = getKey(n);
    if (pairmap.find(key) == pairmap.end()) {
      cout << "missing key '" << key.c_str() <<"'"<<endl;
      cout << n[0] << " " << n[1] << " " << n[2] << " " << n[3] << " " << n[4] << " " << n[5] << endl;
      return;
    }
    bool * b2 = pairmap[key];

    std::vector<std::pair<double, double> > lifespans;
    infile >> s >> n[0];
    assert(strcmp(s,"life") == 0);
    double d[4];
    for (int i=0; i<n[0]; i++) {
      infile >> d[0] >> d[1] >> d[2] >> d[3];
      lifespans.push_back(make_pair(angle(d[0], d[1]), angle(d[2], d[3])));
    }

    bool * b1 = lifespansToAngleRanges(lifespans, num_rotations);

    cout << "lifespans:";
    for (int i=0; i< lifespans.size(); i++)
      cout << "  " << lifespans[i].first << ", " << lifespans[i].second << endl;

    cout << "angle ranges from lifespans: ";
    for (int i=0; i< num_rotations; i++)
      cout << b1[i];
    cout << endl;

    cout << "angle ranges from approximation:  ";
    for (int i=0; i< num_rotations; i++)
      cout << b2[i];
    cout << endl;

    cout << "result: ";
    for (int i=0; i< num_rotations; i++) {
      if (b1[i])
        lifespan_candidates++;

      if (!b1[i] && !b2[i])
        cout << "0"; 
      else if (!b1[i] && b2[i])
        cout << "0"; //victor can eliminate a candidate I can't
      else if (b1[i] && !b2[i]) {
        cout << "0"; //I can eliminate a candidate victor can't (the goal)
        eliminated_candidates++;
      } else if (b1[i] && b2[i])
        cout << "1"; //the candidate can't be eliminated
    }
    cout << endl << "\%:" << ((double)eliminated_candidates)/lifespan_candidates*100 << endl;
  }
}

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "not enough input arguments" << endl;
    return 0;
  }
  char * directory = argv[1];
  std::string dir(directory);
  int num_rotations = 40;

  int samples_per_rotation = 10;
  PTR<Object<Parameter> > tan_half_angle = new InputParameter(tanHalfAngle(num_rotations));
  PTR<Object<Parameter> > sample_tan_half_angle = new InputParameter(sampleTanHalfAngle(num_rotations, samples_per_rotation));
  Parameter::enable();

  if (GENERATE_INPUT_FILE) generateInputFile(directory); //create pairs.txt

  if (GENERATE_ANGLE_RANGES_FILE) { //create out-xx.txt
    PTR<Point> sin_cos_alpha = new SinCosAlpha(tan_half_angle);
    PTR<Point> sample_sin_cos_alpha = new SinCosAlpha(sample_tan_half_angle);
    
    Polyhedron * robot = loadPoly((dir + "/robot").c_str());
    Polyhedron * obstacle = loadPoly((dir + "/obstacle").c_str());
    
    //read/create blockspace
    std::vector<Polyhedron*> blockspaces;
    if (!BLOCKSPACES_FROM_FILE) {
      FreeSpace * fs  = new FreeSpace(robot, obstacle, tan_half_angle, num_rotations, false);
      blockspaces = fs->blockspaces;
    } else {
      char s[50];
      for (int i=0; i<num_rotations; i++) {
        sprintf(s, "sum%02d-out", i);
        Polyhedron * p = loadPoly(s);
        blockspaces.push_back(p);
      }
    }
    for (int i=0; i<blockspaces.size(); i++)
      blockspaces[i]->computeWindingNumbers();

    //read feature pairs
    std::vector<std::pair<PTR<Feature>, PTR<Feature> > > feature_pairs;
    ifstream infile ((dir + "/pairs.txt").c_str());
    int n[6];
    while (infile >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5]) {
      PTR<Feature> fobs, frob;
      if (n[0] == 3) {
        fobs = new Feature(obstacle->vertices[n[1]]->getP(), 
                           obstacle->vertices[n[2]]->getP(), 
                           obstacle->vertices[n[3]]->getP());
        frob = new Feature(robot->vertices[n[5]]->getP());
      } else if (n[0] == 2) {
        fobs = new Feature(obstacle->vertices[n[1]]->getP(), 
                           obstacle->vertices[n[2]]->getP());
        frob = new Feature(robot->vertices[n[4]]->getP(), 
                           robot->vertices[n[5]]->getP());
      } else if (n[0] == 1) {
        fobs = new Feature(obstacle->vertices[n[1]]->getP());
        frob = new Feature(robot->vertices[n[3]]->getP(), 
                           robot->vertices[n[4]]->getP(), 
                           robot->vertices[n[5]]->getP());
      } else {
        cout << "invalid pair" << endl;
        return 0;
      }
      feature_pairs.push_back(make_pair(frob, fobs));
    }

    int start_index, end_index;
    if (argc > 3) {
      start_index = atoi(argv[2]);
      end_index = atoi(argv[3]);
    } else {
      start_index = 0;
      end_index = feature_pairs.size()-1;
    }

    int non_candidates = 0;
    int candidates = 0;
    time_t start,mid,end;
    time (&start);

    ofstream out;
    out.open((dir + "/out-" + std::to_string(num_rotations) + ".txt").c_str());

    //generate angle ranges and output file
    for (int i=start_index; i<=end_index; i++) {
      PTR<Feature> frob = feature_pairs[i].first;
      PTR<Feature> fobs = feature_pairs[i].second;
      bool * b = candidatePairs(fobs, frob, sin_cos_alpha, sample_sin_cos_alpha, blockspaces, samples_per_rotation);
      time(&mid);

      cout << i << ":" << endl;
      for (int j=0; j<num_rotations; j++) {
        cout<<b[j];
        out<<b[j];
        if (b[j]) candidates++; else non_candidates++;
      }
      cout << endl;
      out << endl;
      cout<<" "<< ((double)non_candidates)/(candidates+non_candidates)*100 << "% of candidates eliminated"<<endl;
      double dif = difftime(mid, start);
      cout<<"average "<<dif/(non_candidates+candidates)*num_rotations<<" seconds per feature pair" << endl;
    }

    out.close();

    time (&end);
    double dif = difftime (end,start);
    cout<<"computation time "<<dif<<" seconds"<<endl;
  }

  if (COMPARE_LIFESPANS) {
    testLifeSpans(dir, std::string("out-")+std::to_string(num_rotations)+".txt", num_rotations);
  }
}
