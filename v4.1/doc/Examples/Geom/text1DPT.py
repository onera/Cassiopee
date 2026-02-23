# - text1D (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.text1D("CASSIOPEE")
C.convertPyTree2File(a, 'out.cgns')
