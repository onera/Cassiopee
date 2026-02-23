# - triangulateExteriorFaces (array) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import KCore.test as test

m = C.convertFile2PyTree('boolNG_M1.tp')
m = C.convertArray2NGon(m)

m = XOR.closeCells(m)
C.convertPyTree2File(m, 'out.cgns')
