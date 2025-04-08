# - setInterpData (pyTree)-
# cas structure double wall
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder( (0,0,0), 1, 2, 0, 360, 1, (60, 20, 3) )
b = G.cylinder( (0,0,0), 1, 2, 3, 160, 1, (30, 20, 3) )
C._addBC2Zone(a, 'wall', 'FamilySpecified:SKN1', 'jmin')
C._addBC2Zone(a, 'match', 'BCMatch', 'imin', a, 'imax', trirac=[1,2,3])
C._addBC2Zone(a, 'match', 'BCMatch', 'imax', a, 'imin', trirac=[1,2,3])
C._addBC2Zone(b, 'wall', 'BCWall', 'jmin')
C._addBC2Zone(b, 'wall', 'FamilySpecified:SKN2', 'jmin')
C._addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
C._addBC2Zone(b, 'overlap', 'BCOverlap', 'imax')
t = C.newPyTree(['Base',a,'Base2',b])
C._fillEmptyBCWith(t,'nref','BCFarfield')
C._initVars(t,'Density',1.); C._initVars(t,'centers:G',10.)
C._addState(t[2][1], 'EquationDimension',2)
C._addFamily2Base(t[2][1], 'SKN1', bndType="BCWall")
C._addFamily2Base(t[2][2], 'SKN2', bndType="BCWall")

t1 = X.applyBCOverlaps(t, depth=2)
t1[2][2] = X.setInterpData(t1[2][2],t1[2][1], loc='centers')
test.testT(t1,1)
