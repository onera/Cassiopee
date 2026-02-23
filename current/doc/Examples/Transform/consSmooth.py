# - consSmooth (array) -
import Transform as T
import Converter as C
import Geom as D
a = D.circle((0,0,0), 1., N=20)
a = C.convertArray2Tetra(a)
b = T.consSmooth(a)
C.convertArrays2File([a,b], "out.plt")
