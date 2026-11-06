# - optimizeOverlap (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import KCore.test as test

a1 = G.cylinder((0.,0.,0.), 0.5, 1., 0., 180., 1.,(25,25,10)); a1[0] = 'mesh1'
a2 = G.cylinder((0.,0.,0.), 0.5, 1., 45, 90., 1.,(50,50,10)); a2[0] = 'mesh2'
a1 = C.addBC2Zone(a1,'wall1','BCWall','imin')
a1 = C.addBC2Zone(a1,'wall2','BCWall','imax')
a1 = C.addBC2Zone(a1,'wall3','BCWall','jmin')
a1 = C.fillEmptyBCWith(a1,'nref','BCFarfield')
a2 = C.addBC2Zone(a2,'wall1','BCWall','jmin')
a2 = C.addBC2Zone(a2,'ov1','BCOverlap','imin')
a2 = C.addBC2Zone(a2,'ov2','BCOverlap','imax')
a2 = C.fillEmptyBCWith(a2,'nref','BCFarfield')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a1); t[2][2][2].append(a2)
t[2][1] = C.addState(t[2][1], 'EquationDimension',3)
t[2][2] = C.addState(t[2][2], 'EquationDimension',3)

t = X.applyBCOverlaps(t, depth=1)
t2 = X.optimizeOverlap(t, double_wall=1)
test.testT(t2)

# Same as first test, but now using an intersection dictionnary
interDict = X.getIntersectingDomains(t, method='hybrid')
t2 = X.optimizeOverlap(t, double_wall=1, intersectionsDict=interDict)
test.testT(t2,2)
