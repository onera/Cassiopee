# - rmBCOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = C.rmBCOfType(a, 'BCWall')
t = C.newPyTree(['Base',a])
test.testT(t,1)

# sur un arbre
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imin')
t = C.newPyTree(['Base',a,b])
t = C.rmBCOfType(t, 'BCWall')
test.testT(t,2)

# Sur une liste de zones
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imin')
A = C.rmBCOfType([a,b], 'BCWall')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t,3)

# BC referencees par des familles
a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
a = C.addBC2Zone(a, 'nref', 'FamilySpecified:LOIN', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a]
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')
t[2][1] = C.addFamily2Base(t[2][1], 'LOIN', bndType='BCFarfield')
t = C.rmBCOfType(t, 'BCWall')
test.testT(t,4)

# Rm de familles
a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
a = C.addBC2Zone(a, 'nref', 'FamilySpecified:LOIN', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a]
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')
t[2][1] = C.addFamily2Base(t[2][1], 'LOIN', bndType='BCFarfield')
t = C.rmBCOfType(t, 'FamilySpecified:CARTER')
test.testT(t,5)
#
