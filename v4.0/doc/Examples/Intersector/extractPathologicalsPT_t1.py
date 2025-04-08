# - Extract pathological cells (uncomputable or non-star) (PyTree) -
import Converter.PyTree as C
import Intersector.PyTree as XOR
import KCore.test as test

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

t = XOR.booleanMinus(M1, M2, tol, preserve_right=0, solid_right=0, agg_mode=1)

t = XOR.extractPathologicalCells(t, 2) # ask for 2 level of neighgbors

test.testT(t, 1)
