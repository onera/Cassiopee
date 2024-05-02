# - circle (array) -
import Geom as D
import Converter as C

a = D.circle((0,0,0), 1. , 0., 360.)
C.convertArrays2File(a, "out.plt")
