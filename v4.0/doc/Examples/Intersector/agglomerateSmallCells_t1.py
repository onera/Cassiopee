# - agglomerateSmallCells (array) -
import Generator as G
import Geom as D
import Converter as C
import Intersector as XOR
import KCore.test as test

# octree
s = D.sphere((0,0,0), 1., 100); snear = 0.1
t = G.octree([s],[snear], dfar=5., balancing=1, ratio=2)

print("ngon conversion...")
t = C.convertArray2NGon(t)
print("ngon conformization...")
t = C.conformizeNGon(t); t = G.close(t)
print("ngon close cells...")
t = XOR.closeCells(t)
#t = XOR.reorientExternalFaces(t)

s = C.convertArray2NGon(s)
s = XOR.convertNGON2DToNGON3D(s)
#s = XOR.reorientExternalFaces(s)

x = XOR.diffSurf(t, s, tol=0., preserve_right=1, agg_mode=2) # agg_mode=2 : full mode aggregation

x = XOR.agglomerateSmallCells(x, vmin=0., vratio=0.1)

test.testA(x,1)
