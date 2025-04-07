# - Extract pathological cells (uncomputable or non-star) - (array)

import Converter.PyTree as C
import Intersector.PyTree as XOR
import KCore.test as test

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

t = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=1, agg_mode=1)

t=XOR.computeGrowthRatio(t)
test.testT(t, 1)
