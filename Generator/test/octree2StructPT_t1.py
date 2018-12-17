# - octree2Struct (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

# cas 2D : contours->QUAD sans equilibrage
s = D.circle((0,0,0), 1., N=100); snear = 0.2
res = G.octree([s], [snear], dfar=5., balancing=0)
res1 = G.octree2Struct(res, vmin=5, ext=0)
test.testT(res1,1)
res1 = G.octree2Struct(res, vmin=5, ext=0, optimized=0)
test.testT(res1,2)
#
res1 = G.octree2Struct(res, vmin=5, ext=1)
test.testT(res1,3)
res1 = G.octree2Struct(res, vmin=5, ext=1, optimized=0)
test.testT(res1,4)
#
res1 = G.octree2Struct(res, vmin=5, ext=2)
test.testT(res1,5)
res1 = G.octree2Struct(res, vmin=5, ext=2, optimized=0)
test.testT(res1,6)
#
# cas 2D : contours->QUAD avec equilibrage
res = G.octree([s],[snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=0)
test.testT(res1,7)
res1 = G.octree2Struct(res, vmin=5, ext=0, optimized=0)
test.testT(res1,8)
#
res1 = G.octree2Struct(res, vmin=5, ext=1)
test.testT(res1,9)
res1 = G.octree2Struct(res, vmin=5, ext=1, optimized=0)
test.testT(res1,10)
#
res1 = G.octree2Struct(res, vmin=5, ext=2)
test.testT(res1,11)
res1 = G.octree2Struct(res, vmin=5, ext=2, optimized=0)
test.testT(res1,12)
