# - boolean intersection (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C

M1 = C.convertFile2PyTree('boolNG_M1.tp')
M1 = C.convertArray2NGon(M1)

M2 = C.convertFile2PyTree('boolNG_M2.tp')
M2 = C.convertArray2NGon(M2)

tol = 1.e-12
x = XOR.booleanIntersection(M1, M2, tol, preserve_right=1, solid_right=1)
t = C.newPyTree(['Base',2,x])
C.convertPyTree2File(t, 'boolNGinter11.cgns')

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=0, solid_right=1)
C.convertPyTree2File(x, 'boolNGinter01.cgns')

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=1, solid_right=0)
C.convertPyTree2File(x, 'boolNGinter10.cgns')

x = XOR.booleanIntersection(M1, M2, tol, preserve_right=0, solid_right=0)
C.convertPyTree2File(x, 'boolNGinter00.cgns')
