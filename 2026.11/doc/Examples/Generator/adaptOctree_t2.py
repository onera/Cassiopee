# - adaptOctree (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=5.,balancing=1)
indic = C.node2Center(o)
indic = C.initVars(indic,'indicator',2.)
res = G.adaptOctree(o, indic)
test.testA([res])

# 3D
s = D.sphere((0,0,0), 1., N=20); snear = 0.5
o = G.octree([s], [snear], dfar=5.,balancing=1)
indic = C.node2Center(o)
indic = C.initVars(indic,'indicator',2.)
res = G.adaptOctree(o, indic)
test.testA([res],2)
