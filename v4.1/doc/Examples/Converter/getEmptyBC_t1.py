# - getEmptyBC (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Sur une zone
a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 2))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'imax')
wins = C.getEmptyBC(a,2)
test.testO(wins, 1)

# Sur un arbre
b = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 2))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
wins = C.getEmptyBC(t, 2)
test.testO(wins, 2)

# Sur une liste de zones
wins = C.getEmptyBC([a,b], 2)
test.testO(wins, 3)

# Sur un arbre
b = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 2))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'ov', 'BCOverlap', dim=2)
wins = C.getEmptyBC(t,2)
test.testO(wins,4)
