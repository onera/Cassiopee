# - octree (array) -
# Test avec ratio 3
import Generator as G
import Geom as D
import KCore.test as test

# cas 2D: contours->QUAD
s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5., balancing=0, ratio=3)
test.testA([res],1)
s2 = D.circle((0,0,0), 1., N=100); snear = 0.1
res2 = G.octree([s2],[snear], dfar=5., balancing=1, ratio=3)
test.testA([res2],2)
# cas 3D: surface TRI->HEXA
s = D.sphere((0,0,0), 1., 100); snear = 0.1
res = G.octree([s],[snear], dfar=5., balancing=0, ratio=3)
test.testA([res],3)
s = D.sphere((0,0,0), 1., 100); snear = 0.1
res = G.octree([s],[snear], dfar=5., balancing=1, ratio=3)
test.testA([res],4)
test.writeCoverage(100)
