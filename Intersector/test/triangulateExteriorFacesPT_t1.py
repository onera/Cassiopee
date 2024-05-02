# - triangulateExteriorFaces (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import KCore.test as test

t = C.convertFile2PyTree('boolNG_M1.tp')
C._initVars(t, '{F}={CoordinateX}')
C._initVars(t, '{centers:G}={centers:CoordinateX}')
t = C.convertArray2NGon(t)

t = XOR.triangulateExteriorFaces(t, improve_qual=0)
test.testT(t, 1)
