# - compIndicatorValue(pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Post.PyTree as P

s = D.circle((0,0,0),1.)
snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=1)
res = G.octree2Struct(o, vmin=11,merged=1)
res = G.getVolumeMap(res)
o = P.computeIndicatorValue(o,res,'centers:vol')
t = C.newPyTree(['Base']); t[2][1][2] +=[o]
C.convertPyTree2File(t,"out.cgns")
