# - dist2WallsEikonal (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Converter.Internal as Internal
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test
import numpy

DEPTH = 2
# Bloc cartesien
N = 128; h = 0.1
a = G.cart((0.,0.,0.),(h,h,h),(N,N,1))

# Init wall
sphere = D.sphere((6.4,6.4,0), 1., 100)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base', a])
C._initVars(t,'cellN', 1)
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='node_in')
# Condition aux limites
t = X.setHoleInterpolatedPoints(t,depth=1,loc='nodes')
C._initVars(t,'{flag}=({cellN}>1.)')
t = DTW.distance2WallsEikonal(t,sphere,tc=None,DEPTH=DEPTH,nitmax=10)
test.testT(t,1)

# En centres
# Bloc cartesien
N = 128; h = 0.1
a = G.cart((0.,0.,0.),(h,h,h),(N,N,1))

# Init wall
sphere = D.sphere((6.4,6.4,0), 1., 100)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
C._initVars(t,'centers:cellN', 1)
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='center_in')
# Condition aux limites
t = X.setHoleInterpolatedPoints(t,depth=1,loc='centers')
C._initVars(t,'{centers:flag}=({centers:cellN}>1.)')
t = DTW.distance2WallsEikonal(t,sphere,tc=None,DEPTH=DEPTH,loc='centers',nitmax=10)
test.testT(t,2)
