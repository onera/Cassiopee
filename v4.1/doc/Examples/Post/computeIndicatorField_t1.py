# - compIndicatorField (array) -
import Generator as G
import Converter as C
import Geom as D
import Post as P
import KCore.test as test

s = D.circle((0,0,0), 1., N=100); snear = 0.1
octree = G.octree([s], [snear], dfar=10.,balancing=1)
npts = len(octree[1][0])
indicVal = G.getVolumeMap(octree)

# test sans bodies
indicator, valInf, valSup = P.computeIndicatorField(octree, indicVal)
test.testA([indicator],1)

# test avec bodies
indicator, valInf, valSup = P.computeIndicatorField(
    octree, indicVal, nbTargetPts=2.*npts, bodies=[s])
test.testA([indicator],2)
