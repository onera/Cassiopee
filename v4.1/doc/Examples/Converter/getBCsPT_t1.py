# - getBCs (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(10,10,2))
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
a = C.fillEmptyBCWith(a, 'wall', 'BCWall', dim=2)
(BCs,BCNames,BCTypes) = C.getBCs(a)
test.testO(BCNames, 1)

a = G.cartHexa((0,0,0),(1,1,1),(10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
(BCs,BCNames,BCTypes) = C.getBCs(a)
test.testO(BCNames, 2)

a = G.cartNGon((0,0,0),(1,1,1),(10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
(BCs,BCNames,BCTypes) = C.getBCs(a)
test.testO(BCNames, 3)
