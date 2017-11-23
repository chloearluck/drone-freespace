#include "geometry3d.h"
#include <cstring>
#include <queue>
#include <map>

void findPath(Polyhedron * blockspace, int cell_index, PTR<Point> start, PTR<Point> end, bool startOnSurface, bool endOnSurface, Points &path);
void findPath(Polyhedron * blockspace, PTR<Point> start, PTR<Point> end, Points &path);

