# - setDoublyDefinedBC (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test
ni = 20
dh = 10./ni

a = G.cart((0,0,0),(dh,dh,dh),(ni,ni,ni))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
a = C.addBC2Zone(a, 'overlap_dd', 'BCOverlap', 'kmin',[b],'doubly_defined')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmax')
a = C.addBC2Zone(a, 'overlap3', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'overlap4', 'BCOverlap', 'imax')

t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)

t = C.initVars(t, 'centers:cellN', 1)
t = X.applyBCOverlaps(t)
t = X.setDoublyDefinedBC(t)
test.testT(t)
