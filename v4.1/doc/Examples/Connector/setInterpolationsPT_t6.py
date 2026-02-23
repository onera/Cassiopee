# - setInterpolations (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Transform.PyTree as T
import KCore.test as test

LOCAL = test.getLocal()

NK = 2
DIM = 3
if NK == 2: DIM = 2
a = G.cylinder( (0,0,0), 1, 2, 10, 130, 1, (60, 20, NK)); a[0] = 'cyl1'
b = G.cylinder( (0,0,0), 1, 1.5, 0, 30, 1, (30, 20, NK)); b[0] = 'cyl2'
C._addBC2Zone(a, 'wall', 'BCWall', 'jmin')
C._addBC2Zone(b, 'wall', 'BCWall', 'jmin')
C._fillEmptyBCWith(b,'overlap','BCOverlap',dim=DIM)
t = C.newPyTree(['Base','Base2']); t[2][1][2] = [a]; t[2][2][2] = [b]
C._initVars(t,'Density',1.); C._initVars(t,'centers:cellN',1)
t = X.connectMatchPeriodic(t,rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,120.], dim=DIM)
C._fillEmptyBCWith(t,'nref','BCFarfield',dim=DIM)
C._addState(t[2][1], 'EquationDimension',DIM)
t = X.applyBCOverlaps(t, depth=1)
t = X.setInterpolations(t,double_wall=1,storage='direct',prefixFile=LOCAL+'/chm')
#C.convertPyTree2File(t, "out.cgns")
test.testT(t,1)
# test avec Chimere periodique
NK = 51
a = G.cylinder( (0,0,0), 1, 2, 10, 130, 4., (60, 20, NK)); a[0] = 'cyl1'
b = G.cart((0.4,1.2,-0.3),(0.04,0.04,0.1),(11,11,21))
a = X.connectMatchPeriodic(a,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,120.])
C._fillEmptyBCWith(a, "wall", "BCWall")
C._addBC2Zone(a, 'overlapdd', 'BCOverlap', 'kmin', zoneDonor=[b],rangeDonor='doubly_defined')#dd

C._addBC2Zone(b, 'overlap', 'BCOverlap', 'kmax')
C._fillEmptyBCWith(b, "wall", "BCWall")
for rangel in ['imin','imax','jmin','jmax']:
    C._addBC2Zone(b,'overlapdd','BCOverlap',rangel,zoneDonor=[a],rangeDonor='doubly_defined')

import Converter.elsAProfile as CE
CE._addPeriodicDataInSolverParam(b,rotationCenter=[0.,0.,0.],
                                 rotationAngle=[0.,0.,1.],
                                 NAzimutalSectors=3, isChimera=True)
T._rotate(b,(0,0,0),(0,0,1),60.)

t = C.newPyTree(['Base',a,'Base2',b])
C._initVars(t,'Density',1.); C._initVars(t,'centers:cellN',1)
t = X.applyBCOverlaps(t, depth=1)
t = X.setDoublyDefinedBC(t, depth=1)
t = X.setInterpolations(t,double_wall=1,storage='direct',prefixFile=LOCAL+'/chm')
test.testT(t,2)
