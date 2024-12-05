import Generator as G
import Transform as T
import Converter as C
import Converter as I
import Intersector as XOR
import KCore.test as test
import Post as P

t1 = G.cart((0,0,0), (1,1,1), (10,10,10))
t1 = C.convertArray2NGon(t1); t1 = G.close(t1)
t2 = G.cart((1.,0,0), (1,1,1), (10,10,10))
t2 = C.convertArray2NGon(t2); t2 = G.close(t2)

# test 1 : volume/volume
res = XOR.getOverlappingFaces(t1, t2, RTOL=0.05, amax=0.1)

# create a list of polygon list (t1), one list per zone

t = XOR.agglomerateCellsWithSpecifiedFaces(t1, res[0])

C.convertArrays2File([t], "out.plt")

#test 2 : volume/surface

t2 = P.exteriorFaces(t2)
t2 = XOR.convertNGON2DToNGON3D(t2)

res = XOR.getOverlappingFaces(t1, t2, RTOL=0.05, amax=0.1)

t = XOR.agglomerateCellsWithSpecifiedFaces(t1, res[0])

C.convertArrays2File([t], "out1.plt")