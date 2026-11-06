# - chimeraInfo (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test
a = G.cylinder((0,0,0),1.,3.,360,0,1,(200,30,4)); a[0] = 'cylindre1'
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(a, 'ov1', 'BCOverlap', 'jmax')
b = G.cylinder((4,0,0),1.,3.,360,0,1,(200,30,4)); b[0] = 'cylindre2'
C._addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(b, 'ov1', 'BCOverlap', 'jmax')
t = C.newPyTree(['Corps1', 'Corps2'])
t[2][1][2].append(a); t[2][2][2].append(b)
t = X.connectMatch(t, dim=3)
C._fillEmptyBCWith(t,'nref','BCFarfield', dim=3)
C._addState(t, 'EquationDimension', 3)
C._initVars(t,'F', 0.); C._initVars(t,'centers:G',1.)
t = X.applyBCOverlaps(t, depth=1)
t1 = X.setInterpolations(t, loc='cell', storage='direct')
t1 = X.chimeraInfo(t1,type='interpolated')
t1 = X.chimeraInfo(t1,type='extrapolated')
t1 = X.chimeraInfo(t1,type='orphan')
t1 = X.chimeraInfo(t1,type='cellRatio')
t2 = X.chimeraInfo(t1,type='donorAspect')
test.testT(t2,1)
# in place
X._chimeraInfo(t1,type='donorAspect')
test.testT(t1,1)

t2 = X.setInterpolations(t, loc='cell',storage='inverse')
X._chimeraInfo(t2,type='interpolated')
X._chimeraInfo(t2,type='extrapolated')
X._chimeraInfo(t2,type='orphan')
X._chimeraInfo(t2,type='cellRatio')
X._chimeraInfo(t2,type='donorAspect')
test.testT(t2,2)
