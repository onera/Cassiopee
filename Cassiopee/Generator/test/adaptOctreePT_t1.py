# - adaptOctree (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test
# Quadtree
s = D.circle((0,0,0), 1., N=100); snear=0.1
o = G.octree([s], [snear], dfar=5.,balancing=1)
o = C.initVars(o, 'centers:indicator', 1.)
res = G.adaptOctree(o)
t = C.newPyTree(['OCTREE',2]); t[2][1][2] += [res]
test.testT(t)
# Octree
s = D.sphere((0,0,0), 1., N=100); snear=0.5
o = G.octree([s], [snear], dfar=5., balancing=1)
o = C.initVars(o,'centers:indicator',1.)
res = G.adaptOctree(o)
t = C.newPyTree(['OCTREE']); t[2][1][2] += [res]
test.testT(t,2)
# 9-tree
s = D.circle((0,0,0), 1., N=100); snear=0.1
o = G.octree([s], [snear], dfar=5.,balancing=1,ratio=3)
o = C.initVars(o,'centers:indicator',1.)
res = G.adaptOctree(o,ratio=3)
t = C.newPyTree(['OCTREE',2]); t[2][1][2] += [res]
test.testT(t,3)
# 27-tree
s = D.sphere((0,0,0), 1., N=100); snear=0.5
o = G.octree([s], [snear], dfar=5., balancing=1,ratio=3)
o = C.initVars(o,'centers:indicator',1.)
res = G.adaptOctree(o,ratio=3)
t = C.newPyTree(['OCTREE']); t[2][1][2] += [res]
test.testT(t,4)
