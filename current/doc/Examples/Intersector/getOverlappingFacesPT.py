import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as I
import Intersector.PyTree as XOR
import KCore.test as test
import Post.PyTree as P

t1 = G.cart((0,0,0), (1,1,1), (10,10,10))
t1 = C.convertArray2NGon(t1); t1 = G.close(t1)
t2 = G.cart((1.,0,0), (1,1,1), (10,10,10))
t2 = C.convertArray2NGon(t2); t2 = G.close(t2)

# test 1 : volume/volume
res = XOR.getOverlappingFaces(t1, t2, RTOL=0.05, amax=0.1)

# create a list of polygon list (t1), one list per zone
nb_zones = len(res)
t1zones_pgids = []
for i in range(nb_zones):
    t1zones_pgids.append(res[i][0])

t = XOR.agglomerateCellsWithSpecifiedFaces(t1, t1zones_pgids)

C.convertPyTree2File(t, "out.cgns")

#test 2 : volume/surface

t2 = P.exteriorFaces(t2)
t2 = XOR.convertNGON2DToNGON3D(t2)

res = XOR.getOverlappingFaces(t1, t2, RTOL=0.05, amax=0.1)

t1zones_pgids = []
for i in range(nb_zones):
    t1zones_pgids.append(res[i][0])

t = XOR.agglomerateCellsWithSpecifiedFaces(t1, t1zones_pgids)

C.convertPyTree2File(t, "out1.cgns")