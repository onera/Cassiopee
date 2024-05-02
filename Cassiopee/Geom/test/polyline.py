# - polyline (array) -
import Geom as D
import Converter as C

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)])
C.convertArrays2File(a, "out.plt")
