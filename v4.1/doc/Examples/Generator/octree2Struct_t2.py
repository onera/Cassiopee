# - octree2Struct (array) -
import Generator as G
import Geom as D
import KCore.test as test

# cas 3D: surface->HEXA sans equilibrage
s = D.sphere((0,0,0), 1., N=100); snear = 1.
res = G.octree([s], [snear], dfar=5., balancing=0)
res1 = G.octree2Struct(res, vmin=5, ext=0)
test.testA(res1,1)

res1 = G.octree2Struct(res, vmin=5, ext=0, optimized=0)
test.testA(res1,2)

res1 = G.octree2Struct(res, vmin=5, ext=1)
test.testA(res1,3)
res1 = G.octree2Struct(res, vmin=5, ext=1, optimized=0)
test.testA(res1,4)

res1 = G.octree2Struct(res, vmin=5, ext=2)
test.testA(res1,5)
res1 = G.octree2Struct(res, vmin=5, ext=2, optimized=0)
test.testA(res1,6)

# cas 3D: surface->HEXA avec equilibrage
res = G.octree([s],[snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=0)
test.testA(res1,7)
res1 = G.octree2Struct(res, vmin=5, ext=0, optimized=0)
test.testA(res1,8)

res1 = G.octree2Struct(res, vmin=5, ext=1)
test.testA(res1,9)
res1 = G.octree2Struct(res, vmin=5, ext=1, optimized=0)
test.testA(res1,10)

res1 = G.octree2Struct(res, vmin=5, ext=2)
test.testA(res1,11)
res1 = G.octree2Struct(res, vmin=5, ext=2, optimized=0)
test.testA(res1,12)
