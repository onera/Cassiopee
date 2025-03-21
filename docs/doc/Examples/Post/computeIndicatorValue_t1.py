# - compIndicatorValue(array) -
import Generator as G
import Converter as C
import Geom as D
import Post as P
import KCore.test as test
#2D
s = D.circle((0,0,0),1.)
snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=1)
res = G.octree2Struct(o, vmin=11,merged=1)
vol = G.getVolumeMap(res); res = C.node2Center(res)
res = P.computeIndicatorValue(o,res,vol)
test.testA([res])
# 3D
s = D.sphere((0,0,0),1.)
snear = 1.
o = G.octree([s], [snear], dfar=10., balancing=1)
res = G.octree2Struct(o, vmin=11,merged=1)
vol = G.getVolumeMap(res); res = C.node2Center(res)
res = P.computeIndicatorValue(o,res,vol)
test.testA([res],2)
