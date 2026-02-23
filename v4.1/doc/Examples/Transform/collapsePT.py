# - collapse (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cartTetra((0.,0.,0.),(0.1,0.01,1.),(20,2,1))
b = T.collapse(a)
C.convertPyTree2File(b, "out.cgns")
