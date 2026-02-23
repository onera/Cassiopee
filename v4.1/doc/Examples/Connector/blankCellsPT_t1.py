# - blankCells (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test
import numpy

surf = D.sphere((0,0,0), 0.5, 20)
T._rotate(surf,(0.,0.,0.),(0.,1.,0.), 90.)
a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
C._addBC2Zone(a, 'ov', 'BCOverlap', 'jmin')
t = C.newPyTree(['Cart',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
C._fillEmptyBCWith(t, 'wall', 'BCWall')
C._addVars(t, 'Density')
bodies = [[surf]]
C._initVars(t, 'centers:cellN', 1.)

# Matrice de masquage (arbre d'assemblage)
BM = numpy.array([[1]])
t2 = X.blankCells(t, bodies, BM)
test.testT(t2, 1)

# in place
C._initVars(t2, 'centers:cellN', 1.)
X._blankCells(t2, bodies, BM)
test.testT(t2,1)

# masque inverse
BM = numpy.array([[-1]])
t2 = X.blankCells(t, bodies, BM)
test.testT(t2,2)
