# - naca (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.naca(12.)
C.convertPyTree2File(a, 'out.cgns')
