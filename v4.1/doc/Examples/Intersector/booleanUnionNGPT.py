# - boolean union (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = -0.5e-3

x = XOR.booleanUnion(M1, M2, tol, preserve_right=1, solid_right=1)
t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'boolNGunion11.cgns')

x = XOR.booleanUnion(M1, M2, tol, preserve_right=0, solid_right=1)
t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'boolNGunion01.cgns')

x = XOR.booleanUnion(M1, M2, tol, preserve_right=1, solid_right=0)
t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'boolNGunion10.cgns')

x = XOR.booleanUnion(M1, M2, tol, preserve_right=0, solid_right=0)
t = C.newPyTree(['Base',2]); t[2][1][2].append(x)
C.convertPyTree2File(t, 'boolNGunion00.cgns')
