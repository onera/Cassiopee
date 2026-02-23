# - point (array) -
import Geom as D
import Converter as C

a = D.point((0,0,0))
C.convertArrays2File(a, "out.plt")
