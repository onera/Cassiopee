# - boolean minus (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import KCore.test as test

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

x = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=1)
test.testT(x,1)

x = XOR.booleanMinus(M1, M2, tol, preserve_right=0, solid_right=1)
test.testT(x,2)

#~ x = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=0)
#~ test.testT(x,3)
#~
#~ x = XOR.booleanMinus(M1, M2, tol, preserve_right=0, solid_right=0)
#~ test.testT(x,4)
