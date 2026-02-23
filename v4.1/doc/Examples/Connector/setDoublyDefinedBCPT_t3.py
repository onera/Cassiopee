import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test
import Transform.PyTree as T

NK = 51
a = G.cylinder( (0,0,0), 1, 2, 10, 130, 4., (60, 20, NK)); a[0] = 'cyl1'
b = G.cart((0.4,1.2,-0.3),(0.04,0.04,0.1),(11,11,21))
a = X.connectMatchPeriodic(a,rotationCenter=[0.,0.,0.],rotationAngle=[0.,0.,120.])
C._fillEmptyBCWith(a,"wall","BCWall")
C._addBC2Zone(a,'overlapdd','BCOverlap','kmin',zoneDonor=[b],rangeDonor='doubly_defined')#dd
#
C._addBC2Zone(b,'overlap','BCOverlap','kmax')
C._fillEmptyBCWith(b,"wall","BCWall")
for rangel in ['imin','imax','jmin','jmax']:
    C._addBC2Zone(b,'overlapdd','BCOverlap',rangel,zoneDonor=[a],rangeDonor='doubly_defined')
#
t = C.newPyTree(['Base',a,'Base2',b])
C._initVars(t,'Density',1.); C._initVars(t,'centers:cellN',1)
t = X.setDoublyDefinedBC(t,depth=1)
test.testT(t,1)
#
# 2e cas : la fente est chimere periodique
#
b = T.rotate(b,(0,0,0),(0,0,1),30.)
import Converter.elsAProfile as CE
CE._addPeriodicDataInSolverParam(b,rotationCenter=[0.,0.,0.],
                                 rotationAngle=[0.,0.,1.],
                                 NAzimutalSectors=10, isChimera=True)
#
t = C.newPyTree(['Base',a,'Base2',b])
C._initVars(t,'Density',1.); C._initVars(t,'centers:cellN',1)
t = X.setDoublyDefinedBC(t,depth=1)
test.testT(t,2)
