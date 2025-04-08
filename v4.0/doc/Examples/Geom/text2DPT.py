# - text2D (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.text2D("Cassiopee")
C.convertPyTree2File(a, 'out.cgns')
