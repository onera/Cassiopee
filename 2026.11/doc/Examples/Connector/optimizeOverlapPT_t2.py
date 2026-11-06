# - optimizeOverlap (pyTree) -
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
C._initVars(t,'Density',1.); C._initVars(t,'centers:cellN',1)
t[2][1] = C.addState(t[2][1], 'EquationDimension',2)
t2 = X.optimizeOverlap(t, double_wall=1, priorities=['Base2',0,'Base',1])
test.testT(t2)

# Same as first test, but now using an intersection dictionnary
interDict = X.getIntersectingDomains(t, method='hybrid')
t2 = X.optimizeOverlap(t, double_wall=1,priorities=['Base2',0,'Base',1],
                       intersectionsDict=interDict)
test.testT(t2,2)
