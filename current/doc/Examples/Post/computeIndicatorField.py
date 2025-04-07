# - compIndicatorField (array) -
import Generator as G
import Converter as C
import Geom as D
import Post as P

s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=1)
npts = len(o[1][0])
indicVal = G.getVolumeMap(o)
indicator, valInf, valSup = P.computeIndicatorField(
    o, indicVal, nbTargetPts=2.*npts, bodies=[s])
indicator = C.center2Node(indicator)
o = C.addVars([o, indicator])
C.convertArrays2File(o, "out.plt")
