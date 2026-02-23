# - makeDirect (array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
a = T.reorder(a, (1,2,-3)) # indirect now
a = T.makeDirect(a)
C.convertArrays2File(a, 'out.plt')
