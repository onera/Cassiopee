# - setDoublyDefinedBC (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
b = T.splitNParts(b,2)
C._addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmin', zoneDonor=['FamilySpecified:FENTE'], rangeDonor='doubly_defined')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2]+=b
C._addFamily2Base(t[2][2], 'FENTE')
for z in Internal.getZones(t[2][2]): C._tagWithFamily(z,'FENTE')
C._initVars(t, 'centers:cellN', 1)
t = X.applyBCOverlaps(t)
t = X.setDoublyDefinedBC(t)
C.convertPyTree2File(t, 'out.cgns')
