# - boolean modified solid (array) -
import Intersector as XOR
import Converter as C
import KCore.test as test

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

SOLID = C.convertFile2Arrays('boolNG_M2.tp')
SOLID = C.convertArray2NGon(SOLID[0])

tol = -0.5e-3

x = XOR.booleanModifiedSolid(SOLID, M1, tol, preserve_solid=1)
test.testA([x],1)

x = XOR.booleanModifiedSolid(SOLID, M1, tol, preserve_solid=0)
test.testA([x],2)
