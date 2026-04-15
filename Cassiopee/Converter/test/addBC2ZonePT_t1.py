# - addBC2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# - Structure -
# 3D
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
a = C.addBC2Zone(a, 'nearmatch1','BCNearMatch', 'imax', a, 'imin', [1,2,3])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = C.addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmin', \
                 zoneDonor=[a], rangeDonor='doubly_defined')
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# 2D
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,1))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2])
a = C.addBC2Zone(a, 'nearmatch1', 'BCNearMatch', 'imax', a, 'imin', [1,2])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = C.addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmin', \
                 zoneDonor=[a], rangeDonor='doubly_defined')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
test.testT(t, 2)

# 1D
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,1,1))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1])
a = C.addBC2Zone(a, 'nearmatch1', 'BCNearMatch', 'imax', a, 'imin', [1])
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
test.testT(t, 3)

# Sur un arbre
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((5,0.9,0), (0.01,0.01,1.), (20,20,2)); b[0] = 'cart2'
a = C.initVars(a,'F1',1.); a = C.initVars(a,'centers:G1',2.)
b = C.initVars(b,'F2',3.); b = C.initVars(b,'centers:G2',4.)
t = C.newPyTree(['Base']); t[2][1][2] = t[2][1][2] + [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.addBC2Zone(t, 'wall', 'BCWall', 'imin')
test.testT(t, 4)

# Sur une liste de zones
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((5,0.9,0), (0.01,0.01,1.), (20,20,2)); b[0] = 'cart2'
a = C.initVars(a,'F1',1.); a = C.initVars(a,'centers:G1',2.)
b = C.initVars(b,'F2',3.); b = C.initVars(b,'centers:G2',4.)
A = C.addBC2Zone([a,b], 'wall', 'BCWall', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] = t[2][1][2] + A
test.testT(t,5)

# doubly defined par une famille
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = G.cart((2.5,2.5,-2.5),(0.5,0.5,0.5),(10,10,30)); b[0] = 'fente'
C._addBC2Zone(a, 'overlap1', 'BCOverlap', 'kmin',zoneDonor=['FamilySpecified:FENTE'],rangeDonor='doubly_defined')
C._tagWithFamily(b,'FENTE')
t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)
C._addFamily2Base(t[2][2], 'FENTE')
test.testT(t,6)
