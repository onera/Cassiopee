# - getEmptyBC (pyTree) -
# 2D BE
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cartHexa((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 1))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
wins = C.getEmptyBC(a, dim=2)
test.testO(wins, 1)
a = C.addBC2Zone(a, 'wall1', 'BCWall', faceList=wins)
wins = C.getEmptyBC(a, dim=2)
test.testO(wins, 2)

# Sur une liste de zones
b = G.cartTetra((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
wins = C.getEmptyBC([a, b], dim=2)
test.testO(wins, 3)

# Sur un arbre
t = C.newPyTree(['Base']); t[2][1][2] += [a, b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
wins = C.getEmptyBC(t, dim=2)[0]
test.testO(wins, 4)
zones = Internal.getZones(t)
for i, z in enumerate(zones):
    if wins[i]:
        C._addBC2Zone(z, f'inlet{i+1:d}', 'BCFarfield', faceList=wins[i])
test.testT(t, 41)

# Sur un arbre apres appel a fillEmptyBCWith
b = G.cartTetra((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
t = C.newPyTree(['Base']); t[2][1][2] += [a, b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall', dim=2)
wins = C.getEmptyBC(t, dim=2)
test.testO(wins, 5)
