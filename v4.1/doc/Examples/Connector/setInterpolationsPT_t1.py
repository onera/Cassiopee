# - setInterpolations (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

a = G.cylinder((0,0,0),0.,3.,360,0,1,(200,40,2)); a[0] = 'cylindre1'
b = G.cylinder((0,0,0),1.,2.,360,0,1,(200,10,2)); b[0] = 'cylindre2'
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'nref', 'BCFarfield', 'jmax')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'jmin')
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'jmax')
t = C.newPyTree(['Cyl1','Cyl2']); t[2][1][2].append(a); t[2][2][2].append(b)
t = X.connectMatch(t, dim=2)
t = C.fillEmptyBCWith(t,'nref','BCFarfield',dim=2)
for i in range(1,3): t[2][i] = C.addState(t[2][i], 'EquationDimension', 2)
# depth = 2
t1 = X.applyBCOverlaps(t, depth=2)
t1 = X.setInterpolations(t1, loc='cell',storage='direct')
test.testT(t1,1)
# depth = 1
t2 = X.applyBCOverlaps(t, depth=1)
t2 = X.setInterpolations(t2, loc='face',storage='direct')
t2 = X.setInterpolations(t2, loc='cell',storage='direct')
test.testT(t2,2)
