# - refine (array) -
import Geom as D
import Transform as T
import Converter as C

a = D.line((0,0,0), (1,0,0), N=10)
b = D.line((1,0,0), (2,1,0), N=30)
a = T.join([a,b])
a = D.refine(a, N=30)
C.convertArrays2File(a, 'out.plt')
