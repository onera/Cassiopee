# - setDoublyDefinedBC (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmin',[b],'doubly_defined')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)
t[2][1] = C.addState(t[2][1], 'EquationDimension',3)
t[2][2] = C.addState(t[2][2], 'EquationDimension',3)
t = C.initVars(t, 'Density', 1.); t = C.initVars(t, 'centers:Pressure', 1.)
t2 = X.applyBCOverlaps(t, depth=1)
t2 = X.setDoublyDefinedBC(t2, depth=1)
test.testT(t2,1)
t2 = X.applyBCOverlaps(t, depth=2)
t2 = X.setDoublyDefinedBC(t2, depth=2)
test.testT(t2,2)

# GC doubly defined par une famille
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
C._tagWithFamily(b,'FENTE')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)
C._addFamily2Base(t[2][2], 'FENTE')
C._addBC2Zone(t[2][1][2][0], 'overlap1', 'BCOverlap', 'kmin',zoneDonor=['FamilySpecified:FENTE'],rangeDonor='doubly_defined')
C._initVars(t, 'centers:cellN',1.)
t = X.applyBCOverlaps(t)
t = X.setDoublyDefinedBC(t)
test.testT(t,3)
#
# A la CANNELLE
#
import Converter.Internal as Internal
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
C._addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmin',zoneDonor=['FENTE'],rangeDonor='doubly_defined')
ovdd = Internal.getNodesFromName(a,'overlap*')[0]
Internal.setValue(ovdd,None)
ovdd = Internal.getNodeFromName(ovdd,'UserDefinedData')
Internal.createNode('DonorFamilyName','DataArray_t',value='FENTE',parent=ovdd)
C._tagWithFamily(b,'FENTE')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)
C._addFamily2Base(t[2][2], 'FENTE')
C._initVars(t, 'centers:cellN',1.)
t = X.applyBCOverlaps(t)
t = X.setDoublyDefinedBC(t)
test.testT(t,4)
