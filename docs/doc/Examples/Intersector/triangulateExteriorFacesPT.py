# - triangulateExteriorFaces (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C

t = C.convertFile2PyTree('boolNG_M1.tp')
t = C.convertArray2NGon(t)

t = XOR.triangulateExteriorFaces(t)
C.convertPyTree2File(t, 'out.cgns')
