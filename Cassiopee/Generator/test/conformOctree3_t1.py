# - conformOctree3 (array) -
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

s = D.circle((0,0,0), 1., N=100); snear = 1.
res = G.octree([s], [snear], dfar=5., balancing=1, ratio=3)
res = C.initVars(res,'{F}={x}+2*{y}')
res = G.conformOctree3(res)
test.testA([res],1)
