# - convert2elsAxdt (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.elsAProfile as CE
a = G.cylinder((0,0,0),1.,3.,360,0,1,(200,30,2)); a[0] = 'cylindre1'
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(a, 'ov1', 'BCOverlap', 'jmax')
b = G.cylinder((4,0,0),1.,3.,360,0,1,(200,30,2)); b[0] = 'cylindre2'
C._addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(b, 'ov1', 'BCOverlap', 'jmax')
c1 = G.cart((-5.,-7.5,0), (15./200,15./200,1), (100,200,2))
c2 = G.cart((2.425,-7.5,0), (15./200,15./100,1), (100,100,2))
t = C.newPyTree(['Corps1', 'Corps2', 'Bgrd'])
t[2][1][2].append(a); t[2][2][2].append(b);t[2][3][2]+=[c1,c2]
t = X.connectMatch(t, dim=2)
t = X.connectNearMatch(t, dim=2)
t = X.applyBCOverlaps(t, depth=1)

# blanking
C._initVars(t[2][3],'{centers:cellN}={centers:cellN}*(1.-({centers:CoordinateX}<2)*({centers:CoordinateX}>-2)*({centers:CoordinateY}<2)*({centers:CoordinateY}>-2))')
t = X.setInterpolations(t, loc='cell',storage='inverse')
C._initVars(t,"centers:TurbulentDistance", 1.)
CE._convert2elsAxdt(t)
C.convertPyTree2File(t,"out.cgns")
