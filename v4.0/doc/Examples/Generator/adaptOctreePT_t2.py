# - adaptOctree (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=5.,balancing=1)
o = C.initVars(o,'centers:indicator', 2.)
res = G.adaptOctree(o)
t = C.newPyTree(['OCTREE',2]); t[2][1][2] += [res]
test.testT(t)

s = D.sphere((0,0,0), 1., N=20); snear = 0.5
o = G.octree([s], [snear], dfar=5., balancing=1)
o = C.initVars(o,'centers:indicator',2.)
res = G.adaptOctree(o)
t = C.newPyTree(['OCTREE']); t[2][1][2] += [res]
test.testT(t,2)
