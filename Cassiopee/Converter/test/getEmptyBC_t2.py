# - getEmptyBC (pyTree) -
# 3d
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 21))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'imax')
wins = C.getEmptyBC(a)
test.testO(wins, 1)
# Sur un arbre
b = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 21))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
wins = C.getEmptyBC(t)
test.testO(wins, 2)

# Sur une liste de zones
wins = C.getEmptyBC([a,b])
test.testO(wins, 3)

# Sur un arbre
b = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 21))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',21.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'ov', 'BCOverlap')
wins = C.getEmptyBC(t)
test.testO(wins,4)
