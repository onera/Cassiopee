# - compIndicatorValue(pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Post.PyTree as P
import KCore.test as test
#2D
s = D.circle((0,0,0),1.)
snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=1)
res = G.octree2Struct(o, vmin=11,merged=1)
res = G.getVolumeMap(res)
o = P.computeIndicatorValue(o,res,'centers:vol')
test.testT(o)
# 3D
s = D.sphere((0,0,0),1.)
snear = 1.
o = G.octree([s], [snear], dfar=10., balancing=1)
res = G.octree2Struct(o, vmin=11,merged=1)
res = G.getVolumeMap(res)
o = P.computeIndicatorValue(o,res,'centers:vol')
test.testT(o,2)
