# - fillEmpyBCWith (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

#-----------------
# cas 3D structure
#-----------------
# Avec une seule zone
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
t = C.newPyTree(['Base', a])
test.testT(t,1)

# Avec un arbre
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
b = G.cart((10,0,0),(1,1,1),(20,20,20))
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall')
test.testT(t,2)

# Avec une liste de zones
a = G.cart((0,0,0),(1,1,1),(20,20,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
b = G.cart((10,0,0),(1,1,1),(20,20,2))
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
A = C.fillEmptyBCWith([a,b], 'wall', 'BCWall')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t,3)

# Avec une seule zone
a = G.cart((0,0,0),(1,1,1),(10,10,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
t = C.newPyTree(['Base', a])
test.testT(t, 4)

#-----------------
# cas 2D sturcture
#-----------------
# Avec un arbre
a = G.cart((0,0,0),(1,1,1),(10,10,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
b = G.cart((10,0,0),(1,1,1),(10,10,2))
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall',dim=2)
test.testT(t, 5)

# Avec une liste de zones
a = G.cart((0,0,0),(1,1,1),(20,20,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
b = G.cart((10,0,0),(1,1,1),(20,20,2)); b[0] = 'cart2'
b = C.addBC2Zone(b, 'overlap', 'BCOverlap', 'imin')
A = C.fillEmptyBCWith([a,b], 'wall', 'BCWall',dim=2)
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 6)
