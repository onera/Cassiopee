# - extractNonStar (array) -
import Converter as C
import Intersector as XOR

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])

tol = -0.5e-3

m = XOR.booleanMinus(M1, M2, tol, preserve_right=0, solid_right=0, agg_mode=1)

m = XOR.extractPathologicalCells(m)

C.convertArrays2File(m, "out.plt")
