# - boolean union (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import KCore.test as test

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

M1   = C.fillEmptyBCWith(M1, 'wall','BCWall')
M2   = C.fillEmptyBCWith(M2, 'farfield','BCFarfield')

M1 = C.initVars(M1,'{centers:var1} = 3.')
M2 = C.initVars(M2,'{centers:var2} = 4.')
M1 = C.initVars(M1,'{centers:var}  = 5.')
M2 = C.initVars(M2,'{centers:var}  = 6.')

x = XOR.booleanUnion(M1, M2, tol, preserve_right=1, solid_right=1)
test.testT(x,1)

x = XOR.booleanUnion(M1, M2, tol, preserve_right=0, solid_right=1)
test.testT(x,2)
