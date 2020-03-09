# - boolean difference (array) -
import Intersector as XOR
import Converter as C
import Transform as T

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])
M1 = C.conformizeNGon(M1)
M1 = XOR.closeCells(M1)


M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])
M2 = C.conformizeNGon(M2)
M2 = XOR.closeCells(M2)

tol = -0.5e-3


M = T.join(M1,M2)
M = XOR.selfX(M)

C.convertArrays2File([M], 'out.plt')
