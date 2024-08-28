# - quadrangle (array) -
import Geom as D
import Converter as C

a = D.quadrangle((0,0,0.1), (0.1,0.,0.1), (0.05, 0.08, 0.1), (0.02,0.05,0.1))
C.convertArrays2File(a, "out.plt")
