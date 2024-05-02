# - translate (array) -
import Transform as T
import Generator as G
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))
b = T.translate(a, (-1.,0.,0.))
C.convertArrays2File([a,b], 'out.plt')
