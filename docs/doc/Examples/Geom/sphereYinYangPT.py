# - sphereYinYang (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.sphereYinYang((0,0,0), 1., 50)
C.convertPyTree2File(a, "out.cgns")
