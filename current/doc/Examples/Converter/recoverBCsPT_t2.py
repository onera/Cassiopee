# - recoverBCs (pyTree) -
import Converter.PyTree   as C
import Generator.PyTree   as G
import Post.PyTree        as P
import KCore.test         as test

# STRUCT
a = G.cart((0,0,0),(1,1,1),(10,10,2))
C._addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
C._fillEmptyBCWith(a, 'wall', 'BCWall', dim=2)
C._initVars(a,'{centers:varX}={centers:CoordinateY}+{centers:CoordinateX}')
C._initBCDataSet(a,'{var}=1.')

(BCs,BCNames,BCTypes) = C.getBCs(a)

b = C.convertArray2NGon(a, recoverBC=False)
b = P.selectCells(a,'{centers:varX}>10.')
c = b

C._recoverBCs(b,(BCs,BCNames,BCTypes),removeBC=False)
test.testT(b, 1)

C._recoverBCs(c,(BCs,BCNames,BCTypes),removeBC=True)
test.testT(c, 2)

# HEXA
a = G.cartHexa((0,0,0),(1,1,1),(10,10,2))
C._addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
C._fillEmptyBCWith(a, 'sym', 'BCSymmetryPlane', dim=2)
C._initVars(a,'{centers:varX}={centers:CoordinateY}+{centers:CoordinateX}')
C._initBCDataSet(a,'{var}=1.')

(BCs,BCNames,BCTypes) = C.getBCs(a)
b = C.convertArray2NGon(a, recoverBC=False)
b = P.selectCells(a,'{centers:varX}>10.')
c = b

C._recoverBCs(b,(BCs,BCNames,BCTypes),removeBC=False)
test.testT(b, 3)

C._recoverBCs(c,(BCs,BCNames,BCTypes),removeBC=True)
test.testT(c, 4)

# NGON
a = G.cartNGon((0,0,0),(1,1,1),(10,10,2))
C._addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
C._fillEmptyBCWith(a, 'sym', 'BCSymmetryPlane', dim=2)
C._initVars(a,'{centers:varX}={centers:CoordinateY}+{centers:CoordinateX}')
C._initBCDataSet(a,'{var}=1.')

(BCs,BCNames,BCTypes) = C.getBCs(a)
b = P.selectCells(a,'{centers:varX}>10.')
C._addBC2Zone(b, 'wall', 'BCWall', faceList=[1,2])
c = b

C._recoverBCs(b,(BCs,BCNames,BCTypes),removeBC=False)
test.testT(b, 5)

C._recoverBCs(c,(BCs,BCNames,BCTypes),removeBC=True)
test.testT(c, 6)
