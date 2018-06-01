#include "sumfaces.h"
#include "io.h"

const bool BLOCKSPACES_FROM_FILE = true;

void save(std::vector<PTR<Feature> > & features, int n_points, const char * filename) {
  ofstream out;
  out.open(filename);
  out << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS "<<n_points<<" double"<<endl;
  for (int i=0; i<features.size(); i++)
    for (int j=0; j<features[i]->p.size(); j++)
      out << features[i]->p[j]->getApprox().getX().mid() << " " <<
             features[i]->p[j]->getApprox().getY().mid() << " " <<
             features[i]->p[j]->getApprox().getZ().mid() << endl;
  out << "POLYGONS "<<features.size()<<" "<<features.size() + n_points<<endl;
  int index = 0;
  for (int i=0; i<features.size(); i++) {
      out << features[i]->p.size() << " " << index++ << " " << index++ << " " << index++;
      if (features[i]->p.size() == 4)
        out << " " << index++;
      out << endl;
  }
  out.close();
}

void all_pairs(std::vector<std::pair<PTR<Feature>, PTR<Feature> > > & feature_pairs, Polyhedron * robot, Polyhedron * obstacle) {
  for (int i=0; i<robot->vertices.size(); i++)
    for (int j=0; j<obstacle->faces.size(); j++) {
      PTR<Feature> frob = new Feature(robot->vertices[i]->getP());
      Vertices vs;
      obstacle->faces[j]->boundaryVertices(vs);
      assert(vs.size() == 3);
      PTR<Feature> fobs = new Feature(vs[0]->getP(), vs[1]->getP(), vs[2]->getP());
      feature_pairs.push_back(make_pair(frob, fobs));
    }
  for (int i=0; i<robot->faces.size(); i++)
    for (int j=0; j<obstacle->vertices.size(); j++) {
      Vertices vs;
      robot->faces[i]->boundaryVertices(vs);
      assert(vs.size() == 3);
      PTR<Feature> frob = new Feature(vs[0]->getP(), vs[1]->getP(), vs[2]->getP());
      PTR<Feature> fobs = new Feature(obstacle->vertices[j]->getP());
      feature_pairs.push_back(make_pair(frob, fobs));
    }
  for (int i=0; i<robot->edges.size(); i++)
    for (int j=0; j<obstacle->edges.size(); j++) {
      PTR<Feature> frob = new Feature(robot->edges[i]->getH()->getP(), robot->edges[i]->getT()->getP());
      PTR<Feature> fobs = new Feature(obstacle->edges[j]->getH()->getP(), obstacle->edges[j]->getT()->getP());
      feature_pairs.push_back(make_pair(frob, fobs));
    }
}

double tanHalfAngle (int n) {
  return tan((1.0 + 1.0e-8) * M_PI / n);
}

double sampleTanHalfAngle(int num_rotations, int samples_per_rotation) {
  double sample_theta = (2 * M_PI / num_rotations - 1e-8) / samples_per_rotation;
  return tan(sample_theta / 2);
}

int main (int argc, char *argv[]) {
  int num_rotations = 40;
  int samples_per_rotation = 10;
  PTR<Object<Parameter> > tan_half_angle = new InputParameter(tanHalfAngle(num_rotations));
  PTR<Object<Parameter> > sample_tan_half_angle = new InputParameter(sampleTanHalfAngle(num_rotations, samples_per_rotation));

  Parameter::enable();

  Polyhedron * robot = loadPoly("frustum");
  Polyhedron * obstacle = loadPoly("droneRoom");
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

  std::vector<std::pair<PTR<Feature>, PTR<Feature> > > feature_pairs;
  all_pairs(feature_pairs, robot, obstacle);

  //DEBUG
  std::vector<PTR<Feature> > features;
  int n_points = 0;
  for (int i=0; i<feature_pairs.size(); i++) {
    PTR<Feature> f = feature_pairs[i].first->sum(feature_pairs[i].second);
    n_points += f->p.size();
    features.push_back(f);
  }
  save(features, n_points, "features00.vtk");

  PTR<Point> sin_cos_alpha = new SinCosAlpha(tan_half_angle);
  PTR<Point> sample_sin_cos_alpha = new SinCosAlpha(sample_tan_half_angle);
  for (int i=0; i<feature_pairs.size(); i++) {
    PTR<Feature> frob = feature_pairs[i].first;
    PTR<Feature> fobs = feature_pairs[i].second;
    bool * b = candidatePairs(fobs, frob, sin_cos_alpha, sample_sin_cos_alpha, blockspaces, samples_per_rotation);

    for (int j=0; j<40; j++)
      cout<<b[j];
    cout<<endl;
  }
}
