# - splitConnexity (pyTree) -
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D

a = D.text2D("CASSIOPEE")
B = T.splitConnexity(a)
C.convertPyTree2File(B, 'out.cgns')
