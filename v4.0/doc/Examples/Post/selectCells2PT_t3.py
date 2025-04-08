# - selectCells2 (pyTree) -
# Cas test avec un tag en noeuds
import Converter.PyTree  as C
import Generator.PyTree  as G
import Post.PyTree       as P
import KCore.test as test

def F(x, y, z):
    if (x+2*y+z > 20.): return True
    else: return False

# test sur une zone - tag aux noeuds - sans champ en centre
a = G.cart((0,0,0),(1,1,1),(11,11,11))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = C.initVars(a, 'tag', F, ['CoordinateX','CoordinateY','CoordinateZ'])
b = P.selectCells2(a, 'tag')
test.testT(b, 1)
