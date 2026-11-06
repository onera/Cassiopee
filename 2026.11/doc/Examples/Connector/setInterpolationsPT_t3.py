# - setInterpolations (pyTree)-
# Cas double wall
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder( (0,0,0), 1, 2, 0, 360, 1, (60, 20, 2) )
b = G.cylinder( (0,0,0), 1, 2, 3, 160, 1, (30, 20, 2) )
a = C.addBC2Zone(a, 'wall', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imin', a, 'imax', trirac=[1,2,3])
a = C.addBC2Zone(a, 'match', 'BCMatch', 'imax', a, 'imin', trirac=[1,2,3])
a = C.addBC2Zone(a, 'nref', 'BCFarfield', 'jmax')
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.fillEmptyBCWith(b,'overlap','BCOverlap',dim=2)
t = C.newPyTree(['Base','Base2']); t[2][1][2] = [a]; t[2][2][2] = [b]
t = C.initVars(t,'Density',1.); t = C.initVars(t,'centers:G',10.)
C._addState(t[2][1], 'EquationDimension',2)
# # depth = 2
t1 = X.applyBCOverlaps(t, depth=2)
t1 = X.setInterpolations(t1, loc='cell', double_wall=1,storage='direct')
test.testT(t1,1)
# depth = 1
t2 = X.applyBCOverlaps(t, depth=1)
t2 = X.setInterpolations(t2, loc='face', double_wall=1, storage='direct')
t2 = X.setInterpolations(t2, loc='cell', double_wall=1, storage='direct')
test.testT(t2,2)
