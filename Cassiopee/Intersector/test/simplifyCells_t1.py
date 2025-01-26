# - convexify any concave polygon in the mesh (array) -
import Intersector as XOR
import Converter as C
import KCore.test as test

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])

tol = -0.5e-3

m = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=1, agg_mode=1)
#C.convertArrays2File([m], 'i.plt')

m = XOR.simplifyCells(m, 1)

test.testA([m],1)
