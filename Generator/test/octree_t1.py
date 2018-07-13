# - octree (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

# cas 2D : contours->QUAD
s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5.)
test.testA([res], 1)

s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5., balancing=1)
test.testA([res], 2)

# cas 3D : surface TRI->HEXA
s = D.sphere((0,0,0), 1., 100); snear = 0.1
res = G.octree([s], [snear], dfar=5., balancing=0)
test.testA([res], 3)

s = D.sphere((0,0,0), 1., 100); snear = 0.1
res = G.octree([s], [snear], dfar=5., balancing=1)
test.testA([res], 4)
