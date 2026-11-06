# - boolean modified solid (array) -
import Intersector as XOR
import Converter as C

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

SOLID = C.convertFile2Arrays('boolNG_M2.tp')
SOLID = C.convertArray2NGon(SOLID[0])

tol = -0.5e-3

x = XOR.booleanModifiedSolid(SOLID, M1, tol, preserve_solid=1)
C.convertArrays2File([x], 'boolNGmodifiedsolid1.plt')

x = XOR.booleanModifiedSolid(SOLID, M1, tol, preserve_solid=0)
C.convertArrays2File([x], 'boolNGmodifiedsolid0.plt')
