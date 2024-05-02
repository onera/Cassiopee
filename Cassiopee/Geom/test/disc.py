# - disc (array) -
import Geom as D
import Converter as C

a = D.disc((0,0,0), 1.)
C.convertArrays2File(a, 'out.plt')
