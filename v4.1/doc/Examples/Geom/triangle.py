# - triangle (array) -
import Geom as D
import Converter as C

a = D.triangle((0,0,0), (0.1,0.,0.1), (0.05, 0.08, 0.1))
C.convertArrays2File(a, "out.plt")
