# - selectCells (pyTree) -
# Test avec un tag en centres
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

def F(x, y, z): return x+2*y+z

def F2(x):
    if (x > 15.): return True
    else: return False

# test sur une zone + tag aux centres
a = G.cart((0,0,0),(1,1,1),(11,11,11))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.node2Center(a, 'Density')
a = C.initVars(a, 'centers:tag', F2, ['centers:Density'])
b = P.selectCells2(a, 'centers:tag')
t = C.newPyTree(['Base']); t[2][1][2] += [b]
test.testT(t, 1)

# test sur une base
t = C.newPyTree(['Base']); t[2][1][2] += [a]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
b = P.selectCells2(t[2][1], 'centers:tag')
test.testT(b,2)

# test sur un arbre
t = P.selectCells2(t, 'centers:tag')
test.testT(t,3)
