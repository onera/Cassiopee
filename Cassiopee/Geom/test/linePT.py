# - line (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.line((0,0,0), (1,0,0))
C.convertPyTree2File(a, 'out.cgns')
