# - mergeCart (array) -
# octree 2D
import Generator as G
import Geom as D
import Transform as T
import KCore.test as test

s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5., balancing=1)
res1 = G.octree2Struct(res, vmin=5, ext=0)
res1 = T.mergeCart(res1)
test.testA(res1)
