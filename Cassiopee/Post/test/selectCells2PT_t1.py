# - selectCells2 (pyTree) -
# Cas test avec un tag en noeuds
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

def F(x, y, z):
    if (x+2*y+z > 20.): return True
    else: return False

# test sur une zone
a = G.cart((0,0,0),(1,1,1),(11,11,11))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = C.initVars(a, 'tag', F, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'centers:G', 12.)
b = P.selectCells2(a, 'tag')
test.testT(b,1)

# test sur une base
t = C.newPyTree(['Base']); t[2][1][2] += [a]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
b = P.selectCells2(t[2][1], 'tag')
test.testT(b,2)

# test sur un arbre
t = P.selectCells2(t, 'tag')
test.testT(b,3)
