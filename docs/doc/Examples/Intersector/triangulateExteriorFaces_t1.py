# - triangulateExteriorFaces (array) -
import Intersector as XOR
import Converter as C
import KCore.test as test

m = C.convertFile2Arrays('boolNG_M1.tp')
m = C.convertArray2NGon(m[0])

m = XOR.triangulateExteriorFaces(m)
test.testA([m],1)
