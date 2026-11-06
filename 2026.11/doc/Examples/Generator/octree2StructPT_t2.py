# - octree2Struct (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# cas 3D : contours->HEXA sans equilibrage
s = D.sphere((0,0,0), 1., N=20); snear = 1.
res = G.octree([s],[snear], dfar=5., balancing=0)
t = C.newPyTree(['Base']); t[2][1][2] = [res]
res1 = G.octree2Struct(res, vmin=5, ext=0, merged=0, optimized=0)
test.testT(res1,1)
# ext=1 : overlap
res1 = G.octree2Struct(res, vmin=5, ext=1, merged=0,optimized=0)
test.testT(res1,2)
