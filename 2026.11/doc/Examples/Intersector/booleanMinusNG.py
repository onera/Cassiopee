# - boolean difference (array) -
import Intersector as XOR
import Converter as C

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])

tol = -0.5e-3

x = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=1)
C.convertArrays2File([x], 'boolNGminus11.plt')

x = XOR.booleanMinus(M1, M2, tol, preserve_right=0, solid_right=1)
C.convertArrays2File([x], 'boolNGminus01.plt')

#~ x = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=0)
#~ C.convertArrays2File([x], 'boolNGminus10.plt')
#~
#~ x = XOR.booleanMinus(M1, M2, tol, preserve_right=0, solid_right=0)
#~ C.convertArrays2File([x], 'boolNGminus00.plt')
