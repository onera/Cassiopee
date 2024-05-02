# - dist2WallsEikonal (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Converter.Internal as Internal
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import Generator.PyTree as G
import Transform.PyTree as T
import numpy

DEPTH = 2
# Bloc cartesien
N = 128; h = 0.1
a = G.cart((0.,0.,0.),(h,h,h),(N,N,1))
# Init wall
sphere = D.sphere((6.4,6.4,0), 1., 100)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
#
NIT = 10
for it in range(NIT):
    T._translate(sphere,(0.1*it,0,0))
    C._initVars(t,"cellN",1.)
    t = X.blankCells(t, [[sphere]], numpy.array([[1]]), blankingType='node_in')
    t = X.setHoleInterpolatedPoints(t,depth=1,loc='nodes')
    C._initVars(t,'{flag}=({cellN}>1.)')
    tc = Internal.copyRef(t)
    t = DTW.distance2WallsEikonal(t,sphere,tc=tc,DEPTH=DEPTH,nitmax=10)
C.convertPyTree2File(t, 'out.cgns')
