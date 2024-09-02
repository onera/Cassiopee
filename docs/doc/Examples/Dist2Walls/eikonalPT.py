# - eikonal (pyTree) -
import Dist2Walls.PyTree as Dist2Walls
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Connector.PyTree as X
import numpy

# Bloc cartesien
N = 128; h = 0.1
a = G.cart((0.,0.,0.),(h,h,h),(N,N,1))
C._initVars(a,'cellN', 1.)

# Init wall
sphere = D.sphere((6.4,6.4,0), 1., 100)
sphere = C.convertArray2Tetra(sphere)
sphere = G.close(sphere)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
t = X.blankCellsTri(t, [[sphere]], numpy.array([[1]]), blankingType='node_in')
# Initialise le front
C._initVars(t,'{Phi}=1.e12*({cellN}>0.)')
C._initVars(t,'{speed}=%f'%(1./h))

# Eikonal
t = Dist2Walls.eikonal(t)
C.convertPyTree2File(t, 'out.cgns')
