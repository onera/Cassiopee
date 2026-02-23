# - boolean intersection (array) -
import Intersector as XOR
import Converter as C

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])

tol = 1.e-12

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=1, solid_right=1)
C.convertArrays2File([x], 'boolNGinter11.plt')

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=0, solid_right=1)
C.convertArrays2File([x], 'boolNGinter01.plt')

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=1, solid_right=0)
C.convertArrays2File([x], 'boolNGinter10.plt')

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=0, solid_right=0)
C.convertArrays2File([x], 'boolNGinter00.plt')
