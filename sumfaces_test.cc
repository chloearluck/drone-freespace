#include "sumfaces.h"
#include "io.h"
#include <time.h>
#include <string>

const bool BLOCKSPACES_FROM_FILE = true;
const bool GENERATE_INPUT_FILE = false;

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

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "not enough input arguments" << endl;
    return 0;
  }
  char * directory = argv[1];
  std::string dir(directory);

  if (GENERATE_INPUT_FILE) generateInputFile(directory);

  int num_rotations = 40;
  int samples_per_rotation = 10;
  PTR<Object<Parameter> > tan_half_angle = new InputParameter(tanHalfAngle(num_rotations));
  PTR<Object<Parameter> > sample_tan_half_angle = new InputParameter(sampleTanHalfAngle(num_rotations, samples_per_rotation));

  Parameter::enable();

  Polyhedron * robot = loadPoly((dir + "/robot").c_str());
  Polyhedron * obstacle = loadPoly((dir + "/obstacle").c_str());
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

  PTR<Point> sin_cos_alpha = new SinCosAlpha(tan_half_angle);
  PTR<Point> sample_sin_cos_alpha = new SinCosAlpha(sample_tan_half_angle);
  
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

  for (int i=start_index; i<=end_index; i++) {
    PTR<Feature> frob = feature_pairs[i].first;
    PTR<Feature> fobs = feature_pairs[i].second;
    bool * b = candidatePairs(fobs, frob, sin_cos_alpha, sample_sin_cos_alpha, blockspaces, samples_per_rotation);
    time(&mid);

    for (int j=0; j<40; j++) {
      cout<<b[j];
      out<<b[j];
      if (b[j]) candidates++; else non_candidates++;
    }
    out << endl;
    cout<<" "<< ((double)non_candidates)/(candidates+non_candidates)*100 << "% of candidates eliminated"<<endl;
    double dif = difftime(mid, start);
    cout<<"average "<<dif/(non_candidates+candidates)<<" seconds per feature pair" << endl;
  }

  out.close();

  time (&end);
  double dif = difftime (end,start);
  cout<<"computation time "<<dif<<" seconds"<<endl;
}
