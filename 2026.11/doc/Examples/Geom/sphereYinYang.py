# - sphereYinYang (array) -
import Geom as D
import Converter as C

a = D.sphereYinYang((0,0,0), 1., 50)
C.convertArrays2File(a, "out.plt")
