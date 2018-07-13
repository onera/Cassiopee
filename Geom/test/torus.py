# - torus (array) -
import Geom as D
import Converter as C

a = D.torus((0,0,0), 5., 2.)
C.convertArrays2File(a, "out.plt")
