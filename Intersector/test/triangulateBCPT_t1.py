# - triangulateExteriorFaces (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import KCore.test as test

t = C.convertFile2PyTree('boolNG_M1.tp')

t = C.convertArray2NGon(t)

t = C.fillEmptyBCWith(t, 'wall', 'BCWall', dim=3)

XOR._triangulateBC(t, 'BCWall')

test.testT(t, 1)
