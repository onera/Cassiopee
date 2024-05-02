# - octree (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

# 2D : QUADTREE
s = D.circle((0,0,0), 1., N=100); snear = 0.1
s = C.addVars(s,'ro'); s = C.initVars(s,'centers:cellN',1.)
res = G.octree([s], [snear], dfar=5.,ratio=3)
test.testT(res)

# 3D : octree HEXA
s = D.sphere((0,0,0), 1., 100); snear = 0.1
s = C.addVars(s,'ro'); s = C.initVars(s,'centers:cellN',1.)
res = G.octree([s],[snear], dfar=5.,ratio=3)
test.testT(res,2)

# 2D : QUADTREE
s = D.circle((0,0,0), 1., N=100); snear = 0.1
s = C.addVars(s,'ro'); s = C.initVars(s,'centers:cellN',1.)
res = G.octree([s], [snear], dfar=5., balancing=1,ratio=3)
test.testT(res,3)

# 3D : octree HEXA
s = D.sphere((0,0,0), 1., 100); snear = 0.1
s = C.addVars(s,'ro'); s = C.initVars(s,'centers:cellN',1.)
res = G.octree([s],[snear], dfar=5., balancing=1,ratio=3)
test.testT(res,4)
