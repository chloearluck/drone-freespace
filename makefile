CFLAGS = -g -I. -I/opt/local/include -L/opt/local/lib -std=c++11

COMPILE = g++ $(CFLAGS) -c

LINK = g++ $(CFLAGS)

LFLAGS = -lmpfr -lpthread -Lcplex124/lib -lilocplex -lcplex -lconcert

acp.o:	acp.cc acp.h
	$(COMPILE) acp.cc

expander2.o: expander2.C
	$(COMPILE) -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG \
	-DIL_STD -Icplex124/include expander2.C

create_test_example: create_test_example.o simplify.o expander2.o io.o polyhedron.o poly.o ppoly.o triangulate.o object.o acp.o hull.o
	$(LINK) create_test_example.o simplify.o expander2.o io.o polyhedron.o poly.o ppoly.o triangulate.o object.o acp.o hull.o $(LFLAGS) -o create_test_example

create_test_example.o: create_test_example.cc simplify.h expander2.h io.h polyhedron.h triangulate.h octree.h rbtree.h object.h pv.h acp.h hull.h
	$(COMPILE) create_test_example.cc

freespace.o: freespace.cc freespace.h hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) freespace.cc

freespace_test: freespace_test.o freespace.o hull.o polyhedron.o poly.o ppoly.o triangulate.o io.o object.o acp.o geometry3d.o mink.o simplify.o expander2.o
	$(LINK) freespace_test.o freespace.o hull.o polyhedron.o poly.o ppoly.o triangulate.o io.o object.o acp.o geometry3d.o mink.o simplify.o expander2.o $(LFLAGS) -o freespace_test

freespace_test.o: freespace_test.cc freespace.h hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) freespace_test.cc

freespace-graph.o: freespace-graph.cc freespace-graph.h path3d.h polyhedron.h mink.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) freespace-graph.cc

freespace-graph_test.o: freespace-graph_test.cc freespace-graph.h path3d.h polyhedron.h mink.h io.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) freespace-graph_test.cc

freespace-graph_test: freespace-graph_test.o freespace-graph.o path3d.o polyhedron.o poly.o ppoly.o mink.o triangulate.o io.o object.o acp.o geometry3d.o simplify.o expander2.o
	$(LINK) freespace-graph_test.o freespace-graph.o path3d.o polyhedron.o poly.o ppoly.o mink.o triangulate.o io.o object.o acp.o geometry3d.o simplify.o expander2.o $(LFLAGS) -o freespace-graph_test

geometry3d.o: geometry3d.cc geometry3d.h acp.h pv.h object.h polyhedron.h
	$(COMPILE) geometry3d.cc

hull.o: hull.C hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) hull.C

io.o:	io.C io.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) io.C

lattice: lattice.o mink.o simplify.o poly.o ppoly.o expander2.o io.o polyhedron.o triangulate.o object.o acp.o hull.o
	$(LINK) lattice.o mink.o simplify.o poly.o ppoly.o expander2.o io.o polyhedron.o triangulate.o object.o acp.o hull.o $(LFLAGS) -o lattice

lattice.o: lattice.cc mink.h simplify.h expander2.h io.h polyhedron.h triangulate.h octree.h rbtree.h object.h pv.h acp.h hull.h
	$(COMPILE) lattice.cc

mink.o: mink.C mink.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) mink.C

object.o: object.cc object.h pv.h acp.h
	$(COMPILE) object.cc

path3d.o: path3d.cc path3d.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h 
	$(COMPILE) path3d.cc

path3d_test.o: path3d_test.cc path3d.h polyhedron.h io.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) path3d_test.cc

path3d_test: path3d_test.o path3d.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o
	$(LINK) path3d_test.o path3d.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o $(LFLAGS) -o path3d_test

polyhedron.o: polyhedron.C polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) polyhedron.C

poly.o:	poly.cc poly.h ppoly.h object.h pv.h acp.h
	$(COMPILE) poly.cc

ppoly.o: ppoly.cc ppoly.h acp.h
	$(COMPILE) ppoly.cc

simplify.o: simplify.C simplify.h polyhedron.h poly.h io.h octree.h rbtree.h object.h pv.h acp.h expander2.h
	$(COMPILE) simplify.C

sumfaces.o: sumfaces.cc sumfaces.h polyhedron.h poly.h ppoly.h triangulate.h io.h mink.h freespace.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h hull.h geometry3d.h simplify.h
	$(COMPILE) sumfaces.cc

sumfaces_test.o: sumfaces_test.cc sumfaces.h polyhedron.h poly.h ppoly.h triangulate.h io.h mink.h freespace.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h hull.h geometry3d.h simplify.h
	$(COMPILE) sumfaces_test.cc

sumfaces_test: sumfaces_test.o sumfaces.o polyhedron.o poly.o ppoly.o triangulate.o io.o mink.o expander2.o freespace.o object.o acp.o geometry3d.o hull.o simplify.o
	$(LINK) sumfaces_test.o sumfaces.o polyhedron.o poly.o ppoly.o triangulate.o io.o mink.o expander2.o freespace.o object.o acp.o geometry3d.o hull.o simplify.o  $(LFLAGS) -o sumfaces_test

test_mink:  test_mink.o hull.o mink.o simplify.o expander2.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o
	$(LINK) test_mink.o hull.o mink.o simplify.o expander2.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o $(LFLAGS) -o test_mink

test_mink.o: test_mink.cc mink.h simplify.h expander2.h io.h hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) test_mink.cc

test_new_simplify: test_new_simplify.o simplify.o expander2.o mink.o io.o polyhedron.o triangulate.o object.o acp.o
	$(LINK) test_new_simplify.o simplify.o expander2.o mink.o io.o polyhedron.o triangulate.o object.o acp.o $(LFLAGS) -o test_new_simplify

test_new_simplify.o: test_new_simplify.cc simplify.h expander2.h mink.h io.h polyhedron.h triangulate.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) test_new_simplify.cc

triangulate.o: triangulate.C triangulate.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) triangulate.C 

voxelize: voxelize.o polyhedron.o triangulate.o io.o object.o acp.o
	$(LINK) voxelize.o polyhedron.o triangulate.o io.o object.o acp.o $(LFLAGS) -o voxelize

voxelize.o: voxelize.cc polyhedron.h triangulate.h io.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) voxelize.cc

clean: 
	rm -f *.o *~ hull test_mink freespace_test path3d_test freespace-graph_test *.lp
