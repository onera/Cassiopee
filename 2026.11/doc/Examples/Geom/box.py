# - box (array) -
import Geom as D
import Converter as C

a = D.box((0,0,0), (1,1,1))
C.convertArrays2File(a, 'out.plt')
