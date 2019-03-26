CFLAGS = -g -O3 -std=c++11

COMPILE = g++ $(CFLAGS) -c

LINK = g++ $(CFLAGS)

LFLAGS = -lmpfr -lpthread -Lcplex124/lib -lilocplex -lcplex -lconcert

mink:	main.C simplify.o expander2.o poly.o mink.o polyhedron.o triangulate.o io.o \
	object.o acp.o
	$(LINK) main.C simplify.o expander2.o poly.o mink.o polyhedron.o \
	triangulate.o io.o object.o acp.o $(LFLAGS) -omink

pack:	pack.o mink.o polyhedron.o triangulate.o io.o object.o acp.o
	$(LINK) pack.o mink.o polyhedron.o triangulate.o io.o \
	object.o acp.o $(LFLAGS) -opack

cspace:	csmain.C cspace.o poly.o simplify.o expander2.o mink.o polyhedron.o \
	triangulate.o io.o object.o acp.o
	$(LINK) csmain.C cspace.o poly.o simplify.o expander2.o mink.o polyhedron.o \
	triangulate.o io.o object.o acp.o $(LFLAGS) -ocspace

polyhedron.o: polyhedron.C polyhedron.h octree.h object.h pv.h acp.h
	$(COMPILE) polyhedron.C

io.o:	io.C io.h polyhedron.h octree.h object.h pv.h acp.h
	$(COMPILE) io.C

triangulate.o: triangulate.C triangulate.h polyhedron.h octree.h object.h pv.h acp.h
	$(COMPILE) triangulate.C

mink.o: mink.C mink.h polyhedron.h octree.h object.h pv.h acp.h
	$(COMPILE) mink.C

simplify.o: simplify.C simplify.h poly.h ppoly.h polyhedron.h octree.h object.h \
	pv.h acp.h expander2.h
	$(COMPILE) simplify.C

hull.o: hull.C hull.h polyhedron.h object.h pv.h acp.h
	$(COMPILE) hull.C

pack.o:	pack.C pack.h mink.h polyhedron.h octree.h object.h pv.h acp.h
	$(COMPILE) pack.C

cspace.o: cspace.C cspace.h poly.h ppoly.h mink.h polyhedron.h io.h octree.h \
	object.h pv.h acp.h
	$(COMPILE) cspace.C

object.o: object.cc object.h pv.h acp.h
	$(COMPILE) object.cc

acp.o:	acp.cc acp.h
	$(COMPILE) acp.cc

expander2.o: expander2.C expander2.h
	$(COMPILE)  -Wno-ignored-attributes -m64 -fPIC -fno-strict-aliasing \
	-fexceptions -DNDEBUG -DIL_STD -Icplex124/include expander2.C

round.o: round.C round.h simplify.h polyhedron.h object.h acp.h
	$(COMPILE) round.C

hull:	hull.o io.o triangulate.o polyhedron.o object.o acp.o
	$(LINK) hull.o io.o triangulate.o polyhedron.o object.o acp.o \
	$(LFLAGS) -o hull

delaunay.o: delaunay.C delaunay.h polyhedron.h octree.h object.h pv.h acp.h
	$(COMPILE) delaunay.C

delaunay: delaunay.o simplify.o expander2.o poly.o polyhedron.o triangulate.o \
	io.o object.o acp.o
	$(LINK) delaunay.o simplify.o expander2.o poly.o polyhedron.o \
	triangulate.o io.o object.o acp.o \
	$(LFLAGS) -odelaunay

path3d.o: path3d.cc path3d.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h 
	$(COMPILE) path3d.cc

poly_test: poly_test.o acp.o object.o poly.o
	$(LINK) poly_test.o acp.o object.o poly.o -lmpfr -lpthread -o poly_test

poly_test.o: poly_test.cc poly.h object.h pv.h ppoly.h acp.h
	$(COMPILE) poly_test.cc

poly.o:	poly.cc poly.h ppoly.h object.h pv.h acp.h
	$(COMPILE) poly.cc

freespace.o: freespace.cc freespace.h hull.h polyhedron.h octree.h object.h \
	pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) freespace.cc

freespace_test: freespace_test.cc freespace.o polyhedron.o \
	triangulate.o io.o object.o acp.o mink.o simplify.o expander2.o poly.o
	$(LINK) freespace_test.cc freespace.o polyhedron.o \
	triangulate.o io.o object.o acp.o mink.o simplify.o expander2.o poly.o \
	$(LFLAGS) -o freespace_test

freespace-graph.o: freespace-graph.cc freespace-graph.h path3d.h hull.h polyhedron.h octree.h rbtree.h object.h \
	pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) freespace-graph.cc

freespace-graph_test: freespace-graph_test.cc freespace-graph.o freespace.o path3d.o hull.o polyhedron.o \
	triangulate.o io.o object.o acp.o mink.o simplify.o expander2.o poly.o
	$(LINK) freespace-graph_test.cc freespace-graph.o freespace.o hull.o path3d.o polyhedron.o \
	triangulate.o io.o object.o acp.o mink.o simplify.o expander2.o poly.o \
	$(LFLAGS) -o freespace-graph_test

create_test_example.o: create_test_example.cc hull.h polyhedron.h octree.h rbtree.h object.h \
  pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) create_test_example.cc

create_test_example: create_test_example.o hull.o polyhedron.o \
  triangulate.o io.o object.o acp.o mink.o simplify.o expander2.o poly.o
	$(LINK) create_test_example.o hull.o polyhedron.o \
  triangulate.o io.o object.o acp.o mink.o simplify.o expander2.o poly.o \
	$(LFLAGS) -o create_test_example

lattice: lattice.o mink.o simplify.o expander2.o io.o polyhedron.o poly.o triangulate.o object.o acp.o hull.o
	$(LINK) lattice.o mink.o simplify.o expander2.o io.o polyhedron.o poly.o triangulate.o object.o acp.o hull.o $(LFLAGS) -o lattice

lattice.o: lattice.cc mink.h simplify.h expander2.h io.h polyhedron.h triangulate.h octree.h rbtree.h object.h pv.h acp.h hull.h
	$(COMPILE) lattice.cc

clean: 
	rm -f *.o *~ mink hull delaunay pack cspace poly_test freespace_test *.lp
