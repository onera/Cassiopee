# - cloud (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

x, y, z = [0.0, 0.1, 0.2], [0.0, -0.1, -0.2], [0.0, 0.0, 0.0]
a = D.cloud((x, y, z))
C.convertPyTree2File(a, "out.cgns")
