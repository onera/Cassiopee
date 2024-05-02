# - getOversetInfo (pyTree) -
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
C._addVars(t, 'F'); C._initVars(t,'centers:G',1.)
t = X.applyBCOverlaps(t, depth=1)

tDnr = C.node2ExtCenter(t)
t1 = X.setInterpData(t, tDnr, sameName=1, loc='centers',storage='direct')

C._initVars(tDnr,'{F}={CoordinateX}*{CoordinateY}')
C._initVars(t1,'{centers:F}={centers:CoordinateX}*{centers:CoordinateY}')
X._setInterpTransfers(t1,tDnr,cellNVariable='',variables=['F'])
test.testT(t1,0)
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='interpolated')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='extrapolated')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='orphan')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='cellRatio')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='donorAspect')
test.testT(t1,1)
import sys; sys.exit()

tDnr = X.setInterpData(t, tDnr, sameName=1, loc='centers',storage='inverse')
t1 = X.getOversetInfo(t,tDnr, loc='centers', type='interpolated')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='extrapolated')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='orphan')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='cellRatio')
t1 = X.getOversetInfo(t1,tDnr, loc='centers', type='donorAspect')
test.testT(t1,2)

C._initVars(t, 'cellN',2)
tDnr = C.node2ExtCenter(t)
C._initVars(tDnr, 'cellN',1)
t1 = X.setInterpData(t, tDnr, sameName=1, loc='nodes',storage='direct')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='interpolated')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='extrapolated')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='orphan')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='cellRatio')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='donorAspect')
test.testT(t1,3)

tDnr = X.setInterpData(t, tDnr, sameName=1, loc='nodes',storage='inverse')
t1 = X.getOversetInfo(t,tDnr,loc='nodes',type='interpolated')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='extrapolated')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='orphan')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='cellRatio')
t1 = X.getOversetInfo(t1,tDnr,loc='nodes',type='donorAspect')
test.testT(t1,4)
