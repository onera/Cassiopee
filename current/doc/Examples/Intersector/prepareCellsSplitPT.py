# - convexify any concave polygon in the mesh (array) -
import Intersector.PyTree as XOR
import Converter.PyTree as C

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

m = XOR.booleanMinus(M1, M2, tol, preserve_right=1, solid_right=1, agg_mode=2) #full agg to convexify afterward
#C.convertArrays2File([m], 'i.plt')
m = XOR.prepareCellsSplit(m, PH_set=0, split_policy=2, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold=1.e-2)

C.convertPyTree2File(m, 'out.cgns')
