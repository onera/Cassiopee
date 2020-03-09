# - chimeraTransfer (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

a = G.cylinder((0,0,0),1.,3.,360,0,1,(200,30,3)); a[0] = 'cylindre1'
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(a, 'ov1', 'BCOverlap', 'jmax')
b = G.cylinder((4,0,0),1.,3.,360,0,1,(200,30,3)); b[0] = 'cylindre2'
C._addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(b, 'ov1', 'BCOverlap', 'jmax')
c = G.cart((-5.,-7.5,0), (15./200,15./200,1), (200,200,3))
t = C.newPyTree(['Corps1',a,'Corps2',b,'Bgrd',c])
t = X.connectMatch(t, dim=3)
C._fillEmptyBCWith(t,'nref','BCFarfield', dim=3)
X._applyBCOverlaps(t, depth=2)
C._initVars(t,'{F}={CoordinateX}')
t = C.node2Center(t,['F'])
C._initVars(t,'centers:G', 1.)
t1 = X.setInterpolations(t, loc='cell', storage='direct')
t1 = X.chimeraTransfer(t1, storage='direct', variables=['F','centers:F','centers:G'])
test.testT(t1,1)

t1 = X.setInterpolations(t, storage='inverse')
t2,res = X.chimeraTransfer(t1, storage='inverse', variables=['F','centers:F'], mesh='standard')
test.testT(t2,2)

t2,res = X.chimeraTransfer(t1, storage='inverse', variables=['F','centers:F'], mesh='extended')
test.testT(t2,3)
