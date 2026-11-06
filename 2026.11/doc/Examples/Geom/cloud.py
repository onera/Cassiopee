# - cloud (array) -
import Geom as D
import Converter as C

x, y, z = [0.0, 0.1, 0.2], [0.0, -0.1, -0.2], [0.0, 0.0, 0.0]
a = D.cloud((x, y, z))
C.convertArrays2File(a, "out.plt")
