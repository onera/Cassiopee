# - recoverBCs (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0),(1,1,1),(10,10,2))
a = C.addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'jmin', a, 'jmax', [1,2,3])
a = C.fillEmptyBCWith(a, 'wall', 'BCWall', dim=2)
(BCs,BCNames,BCTypes) = C.getBCs(a)
b = C.convertArray2NGon(a)
C._recoverBCs(b,(BCs,BCNames,BCTypes))
C.convertPyTree2File(b,'out.cgns')
