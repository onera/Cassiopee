# - boolean modified solid (PyTree) -
import Intersector.PyTree as G
import Converter.PyTree as C

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

x = XOR.booleanModifiedSolid(M1, M2, tol, preserve_solid=1)
t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'boolNGmodifiedsolid1.cgns')

x = XOR.booleanModifiedSolid(M1, M2, tol, preserve_solid=0)
t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'boolNGmodifiedsolid0.cgns')
