# - expandLayer (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

# 2D
s = D.circle((0.,0.,0.),1.,N=100)
o = G.octree([s], [0.1], dfar=1.,balancing=1)
o2 = G.expandLayer(o, level=0)
test.testT(o2,1)

# 3D
s = D.sphere((0.,0.,0.),1.,N=100)
o = G.octree([s], [0.1], dfar=1.,balancing=1)
o2 = G.expandLayer(o,level=0)
test.testT(o2,2)
