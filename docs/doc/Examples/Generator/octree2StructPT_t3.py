# - octree2Struct AMR (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

# cas 2D : contours->QUAD AMR
s = D.circle((0,0,0), 1., N=100); snear = 0.2
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=0, AMR=1)
test.testT(res1,1)

# cas 3D : surface->HEXA avec equilibrage
s = D.sphere((0,0,0), 1., N=100); snear = 0.5
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=0, AMR=1)
test.testT(res1,2)

# cas 2D : contours->QUAD AMR
s = D.circle((0,0,0), 1., N=100); snear = 0.2
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=1, AMR=1)
test.testT(res1,3)

# cas 3D : surface->HEXA avec equilibrage
s = D.sphere((0,0,0), 1., N=100); snear = 0.5
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=1, AMR=1)
test.testT(res1,4)
