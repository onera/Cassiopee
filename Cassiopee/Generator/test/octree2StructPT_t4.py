# - octree2Struct liste de vmin (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

# cas 2D : contours->QUAD avec liste de vmin
s = D.circle((0,0,0), 1., N=100); snear = 0.26
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0)
test.testT(res1,11)
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, optimized=0)
test.testT(res1,12)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1)
test.testT(res1,13)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1, optimized=0)
test.testT(res1,14)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2)
test.testT(res1,15)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2, optimized=0)
test.testT(res1,16)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, merged=1)
test.testT(res1,17)
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, optimized=0, merged=1)
test.testT(res1,18)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1, merged=1)
test.testT(res1,19)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1, optimized=0, merged=1)
test.testT(res1,110)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2, merged=1)
test.testT(res1,111)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2, optimized=0, merged=1)
test.testT(res1,112)

# cas 3D : contours->QUAD avec liste de vmin
s = D.sphere((0,0,0), 1., N=100); snear = 0.82
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0)
test.testT(res1,21)
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, optimized=0)
test.testT(res1,22)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1)
test.testT(res1,23)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1, optimized=0)
test.testT(res1,24)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2)
test.testT(res1,25)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2, optimized=0)
test.testT(res1,26)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, merged=1)
test.testT(res1,27)
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=0, optimized=0, merged=1)
test.testT(res1,28)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1, merged=1)
test.testT(res1,29)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=1, optimized=0, merged=1)
test.testT(res1,210)

res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2, merged=1)
test.testT(res1,211)
res1 = G.octree2Struct(res, vmin=[5,7,10], ext=2, optimized=0, merged=1)
test.testT(res1,212)
