# - line (array) -
import Geom as D
import Converter as C

a = D.line((0,0,0), (1,0,0))
C.convertArrays2File(a, 'out.plt')
