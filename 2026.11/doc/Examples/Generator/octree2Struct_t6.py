# - octree2Struct (array) -
import Generator as G
import Geom as D
import KCore.test as test

# cas 3D : contours->QUAD avec liste de vmin
s = D.sphere((0,0,0), 1., N=100); snear = 0.5
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0)
test.testA(res1,1)
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, optimized=0)
test.testA(res1,2)
#
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1)
test.testA(res1,3)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1,optimized=0)
test.testA(res1,4)
#
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0,merged=1)
test.testA(res1,7)
# Cas equilibre
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, optimized=0,merged=1)
test.testA(res1,8)
