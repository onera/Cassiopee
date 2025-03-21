# - optimizeOverlap (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import KCore.test as test

H = 1.; NK = 4*int(H)+1
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., H, (11,11,NK))
t = C.newPyTree(['Base','Base2']); t[2][1][2] += [a]
t = X.connectMatchPeriodic(t,rotationCenter=[0.,0.,0.], translation=[0.,0.,H])

H2 = 4*H; NK = 2*int(H2)+1
b = G.cylinder((0.,0.,-H), 0.1, 1., -30., 120., 3*H, (11,11,NK))
t[2][2][2] += [b]

# Same as first test, but now using an intersection dictionnary
interDict = X.getIntersectingDomains(t, method='hybrid')
t1 = X.optimizeOverlap(t, intersectionsDict=interDict)
test.testT(t1,4)
# First test:
t = X.optimizeOverlap(t)
test.testT(t,1)
# Periodicite par rotation
H = 1.; NK = 4*int(H)+1
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., H, (11,11,NK))
t = C.newPyTree(['Base','Base2']); t[2][1][2] += [a]
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           rotationAngle=[0.,0.,90.])
H2 = 4*H; NK = 2*int(H2)+1
b = G.cylinder((0.,0.,-H), 0.1, 1., -30., 120., 3*H, (11,11,NK))
t[2][2][2] += [b]
# Same as second test, but now using an intersection dictionnary
InterDict = X.getIntersectingDomains(t, method='hybrid')
t2 = X.optimizeOverlap(t, intersectionsDict=interDict)
test.testT(t2,5)
# Second test:
t = X.optimizeOverlap(t)
test.testT(t,2)


# cas double wall
a = G.cylinder( (0,0,0), 1, 2, 0, 360, 1, (60, 20, 2) )
b = G.cylinder( (0,0,0), 1, 2, 3, 160, 1, (30, 20, 2) )
C._addBC2Zone(a, 'wall', 'BCWall', 'jmin')
C._addBC2Zone(a, 'match', 'BCMatch', 'imin', a, 'imax', trirac=[1,2,3])
C._addBC2Zone(a, 'match', 'BCMatch', 'imax', a, 'imin', trirac=[1,2,3])
C._addBC2Zone(a, 'nref', 'BCFarfield', 'jmax')
C._addBC2Zone(b, 'wall', 'BCWall', 'jmin')
C._fillEmptyBCWith(b,'overlap','BCOverlap',dim=2)
t = C.newPyTree(['Base',a,'Base2',b])
C._initVars(t,'Density',1.); C._initVars(t,'centers:cellN',1)
C._addState(t[2][1], 'EquationDimension',2)
t = X.connectMatchPeriodic(t,translation=(0,0,1))
# Same as third test, but now using an Intersection Dictionnary
InterDict = X.getIntersectingDomains(t, method='hybrid')
t3 = X.optimizeOverlap(t, double_wall=1, priorities=['Base2',0,'Base',1],
                       intersectionsDict=interDict)
test.testT(t3,6)
# Third test:
t = X.optimizeOverlap(t, double_wall=1, priorities=['Base2',0,'Base',1])
test.testT(t,3)
