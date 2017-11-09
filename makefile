CFLAGS = -g -I/opt/local/include -L/opt/local/lib

COMPILE = g++ $(CFLAGS) -c

LINK = g++ $(CFLAGS)

LFLAGS = -lmpfr -lpthread -Lcplex124/lib -lilocplex -lcplex -lconcert

acp.o:	acp.cc acp.h
	$(COMPILE) acp.cc

expander2.o: expander2.C
	$(COMPILE) -std=c++11 -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG \
	-DIL_STD -Icplex124/include expander2.C

freespace.o: freespace.cc freespace.h hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h simplify.h
	$(COMPILE) freespace.cc

freespace_test: freespace_test.o freespace.o hull.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o mink.o simplify.o expander2.o
	$(LINK) freespace_test.o freespace.o hull.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o mink.o simplify.o expander2.o $(LFLAGS) -o freespace_test

freespace_test.o: freespace_test.cc freespace.h hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) freespace_test.cc

freespace-graph.o: freespace-graph.cc freespace-graph.h path3d.h polyhedron.h mink.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h 
	$(COMPILE) freespace-graph.cc

freespace-graph_test.o: freespace-graph_test.cc freespace-graph.h path3d.h polyhedron.h mink.h io.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) freespace-graph_test.cc

freespace-graph_test: freespace-graph_test.o freespace-graph.o path3d.o polyhedron.o mink.o triangulate.o io.o object.o acp.o geometry3d.o
	$(LINK) freespace-graph_test.o freespace-graph.o path3d.o polyhedron.o mink.o triangulate.o io.o object.o acp.o geometry3d.o $(LFLAGS) -o freespace-graph_test

geometry3d.o: geometry3d.cc geometry3d.h acp.h pv.h object.h polyhedron.h
	$(COMPILE) geometry3d.cc

hull.o: hull.C hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) hull.C

io.o:	io.C io.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) io.C

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

simplify.o: simplify.C simplify.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h expander2.h
	$(COMPILE) simplify.C

test_union:  test_union.o hull.o mink.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o
	$(LINK) test_union.o hull.o mink.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o -lmpfr -o test_union

test_union.o: test_union.cc mink.h hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) test_union.cc

triangulate.o: triangulate.C triangulate.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) triangulate.C

clean: 
	rm -f *.o *~ hull test_union freespace_test path3d_test freespace-graph_test
