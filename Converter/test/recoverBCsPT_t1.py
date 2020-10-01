# - recoverBCs (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(10,10,2))
C._addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
C._fillEmptyBCWith(a, 'wall', 'BCWall', dim=2)
(BCs,BCNames,BCTypes) = C.getBCs(a)
b = C.convertArray2NGon(a)
C._recoverBCs(b,(BCs,BCNames,BCTypes))
test.testT(b, 1)

a = G.cartHexa((0,0,0),(1,1,1),(10,10,2))
C._addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
(BCs,BCNames,BCTypes) = C.getBCs(a)
b = C.convertArray2NGon(a)
C._recoverBCs(b,(BCs,BCNames,BCTypes))
test.testT(b, 2)

a = G.cartNGon((0,0,0),(1,1,1),(10,10,2))
C._addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
(BCs,BCNames,BCTypes) = C.getBCs(a)
C._rmBCOfType(a,'BCWall')
C._recoverBCs(a,(BCs,BCNames,BCTypes))
test.testT(a, 3)


# on ne detruit pas les matches ou BC existantes
a = G.cart((0,0,0),(1,1,1),(20,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
C._addBC2Zone(a, 'overlap', 'BCOverlap', 'imax')
(BCs,BCNames,BCTypes) = C.getBCs(a)
b = C.convertArray2NGon(a,recoverBC=False)
import Transform.PyTree as T
import Connector.PyTree as X
b = T.splitNParts(b,2)
b = X.connectMatch(b)
C._recoverBCs(b,(BCs,BCNames,BCTypes), removeBC=False)
test.testT(b,4)
