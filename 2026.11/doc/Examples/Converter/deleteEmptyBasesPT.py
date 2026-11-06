# - deleteEmptyBases (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (3,3,3))
b = G.cart((2,0,0), (1,1,1), (3,3,3))

t = C.newPyTree(['Base1', a, 'Base2', 'Base3', b])
t = C.deleteEmptyBases(t)
C.convertPyTree2File(t, 'out.cgns')
