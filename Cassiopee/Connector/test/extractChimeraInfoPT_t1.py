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
C._initVars(t,'F',0.); C._initVars(t,'centers:G',1.)
t = X.applyBCOverlaps(t, depth=1)
t1 = X.setInterpolations(t, loc='cell', storage='direct')
X._chimeraInfo(t1,type='interpolated')
X._chimeraInfo(t1,type='extrapolated')
X._chimeraInfo(t1,type='orphan')
interpPts = X.extractChimeraInfo(t1,type='interpolated',loc='centers')
test.testT(interpPts,1)
extrapPts = X.extractChimeraInfo(t1,type='extrapolated',loc='centers')
test.testT(extrapPts,2)
cfExtrapPts = X.extractChimeraInfo(t1,type='cf>1.5',loc='centers')
test.testO(cfExtrapPts,3)
orphanPts = X.extractChimeraInfo(t1,type='orphan',loc='centers')
test.testT(orphanPts,4)
