# - point (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.point((0,0,0))
C.convertPyTree2File(a, "out.cgns")
