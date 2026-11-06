# - addkplane (array) -
import Geom as D
import Transform as T
import Converter as C

a = D.naca(12., 501)
a = T.addkplane(a)
C.convertArrays2File(a, "out.plt")
