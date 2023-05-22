# - blankClosestTargetCells (pyTree) -
import Converter.PyTree as C
import Dist2Walls.PyTree as DTW
import Generator.PyTree as G
import Geom.PyTree as D
import Connector.PyTree as X
import Connector.ToolboxIBM as TIBM

N = 11; N2 = 2*N-1
h = 1./(N-1)
a = G.cart((-1.,-1.,-1.),(h,h,h), (N2,N2,2))
zmean = C.getMeanValue(a,'CoordinateZ')
surf = D.sphere((0,0,zmean), 0.5, 20)
t = C.newPyTree(['Cart',a])
C.convertPyTree2File(surf,"in.cgns")
bodies = [[surf]]
# Matrice de masquage (arbre d'assemblage)
import numpy
BM = numpy.array([ [ 1] ] )
C._initVars(t, 'centers:cellN', 1.)
t = X.blankCells(t, bodies, BM, blankingType='center_in', delta=0.)
DTW._distance2Walls(t,surf,loc='centers',type='ortho')
C._initVars(t,"{centers:F}={centers:cellN}")
TIBM._blankClosestTargetCells(t)
C.convertPyTree2File(t, 'out.cgns')
