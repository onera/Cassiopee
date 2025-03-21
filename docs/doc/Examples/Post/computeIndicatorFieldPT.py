# - compIndicatorField (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Post.PyTree as P

#----------------------
# indicateur en centres
#----------------------
s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=10., balancing=1)
npts = o[1][0][0]
o = G.getVolumeMap(o)
o, valInf, valSup = P.computeIndicatorField(o, 'centers:vol',
                                            nbTargetPts=2*npts, bodies=[s])
C.convertPyTree2File(o, 'out.cgns')
