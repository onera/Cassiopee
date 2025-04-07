# - sphere (array) -
import Geom as D
import Converter as C

a = D.sphere((0,0,0), 1., 20)
C.convertArrays2File(a, "out.plt")
