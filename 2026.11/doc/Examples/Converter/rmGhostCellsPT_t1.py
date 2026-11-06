# - rmGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

a = G.cart((1,1,1), (1.,1.,1.), (4,2,3)); a[0]='cart1'
b = G.cart((4,1,1), (0.5,1.,1.), (4,2,3)); b[0]='cart2'
c = G.cart((1,1,-3), (1.,1.,0.5), (4,2,9)); c[0]='cart3'

a = T.reorder(a, (-2,1,3))
b = T.reorder(b, (1,2,3))
c = T.reorder(c, (3,1,2))

# Physical BC (here BCWall)
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')

# Chimere
#-----------------------------
# initialisation cellNatureField
C._initVars(a, 'centers:cellN', 1.)
C._initVars(b, 'centers:cellN', 0.)
#-----------------------------
# initialisation density
C._initVars(a, '{Density}=3.*{CoordinateX}+2.*{CoordinateY}')
C._initVars(a, 'StagnationPressure', 0.5)
C._initVars(c, 'StagnationPressure', 1.)
t = C.newPyTree(['Base',a,b,c])
# Matching BC
t = X.connectMatch(t)
a = t[2][1][2][0]

ag = Internal.addGhostCells(t, a, 2, adaptBCs=0)
t[2][1][2].append(ag)

ag = Internal.rmGhostCells(t, ag, 1, adaptBCs=0)
t[2][1][2].append(ag)
ag = Internal.rmGhostCells(t, ag, 1, adaptBCs=0)
t[2][1][2].append(ag)

test.testT(t, 1)
