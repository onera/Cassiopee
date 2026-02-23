# - conformOctree3 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

s = D.circle((0,0,0), 1., N=100); snear = 1.
res = G.octree([s], [snear], dfar=5., balancing=1, ratio=3)
res = C.initVars(res,'{F}={CoordinateX}+2*{CoordinateY}')
res = C.node2Center(res,['F'])
res = C.initVars(res,'{centers:G}=3*{centers:F}')
res = C.rmVars(res,['centers:F'])
res = G.conformOctree3(res)
t = C.newPyTree(['Octree']); t[2][1][2] += [res]
test.testT(t,1)
