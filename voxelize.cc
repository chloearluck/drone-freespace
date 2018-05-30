#include "polyhedron.h"
#include "pv.h"
#include "object.h"
#include "io.h"

Polyhedron * loadPoly(const char * filename) {
  Polyhedron * poly;
  ifstream infile (filename);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}


int main (int argc, char *argv[]) {
  std::vector<Polyhedron*> spaces;
  char s[40];
  for (int i=0; i<= 40; i++) {
    sprintf(s, "two_room_output/sum%02d-out.vtk", i);
    // cout<<s<<endl;
    spaces.push_back(loadPoly(s));
    spaces[i]->computeWindingNumbers();
  }

  double center_x = 0.0;
  double center_y = 0.0;
  double center_z = 0.0;

  double r = 7.0;
  int tiles = 16;
  int n = tiles*tiles; //number of samples per dimension
  float incr = 2.0 * r / n; 

  int tmp = 0;
  
  for (int z = 0; z < n; z++) {
    //zth tile
    float z_val = (z - (n-1)/2.0)*incr + center_z;
    for (int y = 0; y < n; y++) {
      float y_val = (y - (n-1)/2.0)*incr + center_y;
      for (int x = 0; x < n; x++) {
        float x_val = (x - (n-1)/2.0)*incr + center_x;
        PTR<Point> p = new InputPoint(x_val, y_val, z_val);
        //to do: find every space which contains p
        cout<<"{";
        for (int i=0; i<40; i++) {
          if (spaces[i]->contains(p))
            cout<<i<<" ";
        }
        cout<<"}"<<endl;


        //for now, just output the coordinates to see if we're doing the correct order
        // cout<<x_val<<" "<<y_val<<" "<<z_val<<endl;
        // tmp++; if (tmp == 500) return 0; //stop after 100 lines
      }


    }
  }

  return 0;

}

