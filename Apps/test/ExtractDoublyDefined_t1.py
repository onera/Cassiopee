# - ExtractDoublyDefined -
import Apps.Chimera.ExtractDoublyDefined as App

import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

# Le gros bloc
a = G.cart((0,0,0),(1,1,1),(10,10,10))
C._fillEmptyBCWith(a, 'wall', 'BCWall')
C._addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmin', zoneDonor=['FamilySpecified:FENTE'], rangeDonor='doubly_defined')

# La fente
b = G.cart((2.,2.,-2.),(0.5,0.5,0.5),(11,11,12)); b[0] = 'fente'
C._addBC2Zone(b, 'wall', 'BCWall', 'kmin')
C._addBC2Zone(b, 'wall', 'BCWall', 'imin')
C._addBC2Zone(b, 'wall', 'BCWall', 'imax')
C._addBC2Zone(b, 'wall', 'BCWall', 'jmin')
C._addBC2Zone(b, 'wall', 'BCWall', 'jmax')
C._addBC2Zone(b, 'overlap', 'BCOverlap', 'kmax')
C._addBC2Zone(b, 'overlap1', 'BCOverlap', 'imin', zoneDonor=['FamilySpecified:BLOC'], rangeDonor='doubly_defined')
C._addBC2Zone(b, 'overlap1', 'BCOverlap', 'imax', zoneDonor=['FamilySpecified:BLOC'], rangeDonor='doubly_defined')
C._addBC2Zone(b, 'overlap1', 'BCOverlap', 'jmin', zoneDonor=['FamilySpecified:BLOC'], rangeDonor='doubly_defined')
C._addBC2Zone(b, 'overlap1', 'BCOverlap', 'jmax', zoneDonor=['FamilySpecified:BLOC'], rangeDonor='doubly_defined')

t = C.newPyTree(['Base1',a,'Base2',b])

C._addFamily2Base(t[2][1], 'BLOC')
C._addFamily2Base(t[2][2], 'FENTE')
for z in Internal.getZones(t[2][2]): C._tagWithFamily(z,'FENTE')
for z in Internal.getZones(t[2][1]): C._tagWithFamily(z,'BLOC')

myApp = App.ExtractDoublyDefined()
myApp.set(dataIn=t)
myApp.run()
