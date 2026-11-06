# - boolean intersection (array) -
import Intersector as XOR
import Converter as C
import KCore.test as test

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])

tol = -0.5e-3

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=1, solid_right=1)
test.testA([x],1)

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=0, solid_right=1)
test.testA([x],2)

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=1, solid_right=0)
test.testA([x],3)

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=0, solid_right=0)
test.testA([x],4)
