# - convexifyFaces (array) -
# convexify any concave polygon in the mesh
import Intersector as XOR
import Converter as C

M1 = C.convertFile2Arrays('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1[0])

M2 = C.convertFile2Arrays('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2[0])

tol = -0.5e-3

m = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=1, agg_mode=2) #full agg to convexify afterward
#C.convertArrays2File([m], 'i.plt')
m = XOR.convexifyFaces(m)

C.convertArrays2File([m], 'out.plt')
