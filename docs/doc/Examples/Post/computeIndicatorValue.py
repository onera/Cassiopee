# - compIndicatorValue(array) -
import Generator as G
import Converter as C
import Geom as D
import Post as P

s = D.circle((0,0,0),1.)
snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=0)
res = G.octree2Struct(o, vmin=11,merged=0)
vol = G.getVolumeMap(res); res = C.node2Center(res)
val = P.computeIndicatorValue(o,res,vol)
o = C.addVars([o,val])
C.convertArrays2File([o], "out.plt")
