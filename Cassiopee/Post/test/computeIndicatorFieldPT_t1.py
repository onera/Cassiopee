# - compIndicatorField (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Post.PyTree as P
import KCore.test as test

#----------------------
# indicateur en centres
#----------------------
s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=1)
o = G.getVolumeMap(o)
# test sans bodies
o, valInf, valSup = P.computeIndicatorField(o, 'centers:vol')
test.testT(o, 1)

# test avec bodies
indicator, valInf, valSup = P.computeIndicatorField(
    o, 'centers:vol', nbTargetPts=5000, bodies=[s])
test.testT(o,2)

# indicateur en noeuds
def initIndic(v,w): return 3*v + 4*v*v

s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=10.,balancing=1)
o = C.initVars(o,'Density',initIndic,['CoordinateX','CoordinateY'])
# test sans bodies
o, valInf, valSup = P.computeIndicatorField(o, 'Density')
test.testT(o,3)
# test avec bodies
indicator, valInf, valSup = P.computeIndicatorField(
    o, 'Density', nbTargetPts=5000, bodies=[s])
test.testT(o,4)
