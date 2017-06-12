CFLAGS = -g -I/opt/local/include -L/opt/local/lib

COMPILE = g++ $(CFLAGS) -c

LINK = g++ $(CFLAGS)

LFLAGS = -lmpfr -lpthread -Lcplex124/lib -lilocplex -lcplex -lconcert

# mink:	main.C simplify.o mink.o polyhedron.o triangulate.o io.o object.o acp.o expander.o
# 	$(LINK) main.C simplify.o mink.o polyhedron.o triangulate.o io.o \
# 	object.o acp.o expander.o $(LFLAGS) -omink

# pack:	pack.o mink.o polyhedron.o triangulate.o io.o object.o acp.o
# 	$(LINK) pack.o mink.o polyhedron.o triangulate.o io.o \
# 	object.o acp.o $(LFLAGS) -opack

# cspace:	cspace.o mink.o polyhedron.o triangulate.o io.o object.o acp.o
# 	$(LINK) cspace.o mink.o polyhedron.o triangulate.o io.o object.o \
#         acp.o $(LFLAGS) -ocspace

# hull:   hull.o polyhedron.o triangulate.o io.o object.o acp.o
# 	$(LINK) hull.o polyhedron.o triangulate.o io.o object.o acp.o -lmpfr -ohull

geometry3d.o: geometry3d.cc geometry3d.h acp.h pv.h object.h polyhedron.h
	$(COMPILE) geometry3d.cc

polyhedron.o: polyhedron.C polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) polyhedron.C

io.o:	io.C io.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) io.C

triangulate.o: triangulate.C triangulate.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) triangulate.C

hull.o: hull.C hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
	$(COMPILE) hull.C

outer_approx:  outer_approx.o hull.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o
	$(LINK) outer_approx.o hull.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o -lmpfr -o outer_approx

outer_approx.o: outer_approx.c hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) outer_approx.c

test_union:  test_union.o hull.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o
	$(LINK) test_union.o hull.o polyhedron.o triangulate.o io.o object.o acp.o geometry3d.o -lmpfr -o test_union

test_union.o: test_union.c hull.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h geometry3d.h
	$(COMPILE) test_union.c

# mink.o: mink.C mink.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
# 	$(COMPILE) mink.C

# simplify.o: simplify.C simplify.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h \
# 	expander.h
# 	$(COMPILE) simplify.C

# pack.o:	pack.C pack.h mink.h polyhedron.h octree.h rbtree.h object.h pv.h acp.h
# 	$(COMPILE) pack.C

# cspace.o: cspace.C cspace.h mink.h polyhedron.h io.h octree.h rbtree.h object.h pv.h acp.h
# 	$(COMPILE) cspace.C

object.o: object.cc object.h pv.h acp.h
	$(COMPILE) object.cc

acp.o:	acp.cc acp.h
	$(COMPILE) acp.cc

# expander.o: expander.C
# 	$(COMPILE) -std=c++11 -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG \
# 	-DIL_STD -Icplex124/include expander.C

clean: 
	rm -f *.o *~ hull test_union outer_approx
