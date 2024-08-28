# - cone (array) -
import Geom as D
import Converter as C

a = D.cone((0,0,0), 1. , 0.5, 1.)
C.convertArrays2File(a, "out.plt")
