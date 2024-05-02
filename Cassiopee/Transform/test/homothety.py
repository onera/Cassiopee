# - homothety (array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))
b = T.homothety(a, (0.,0.,0.), 2.)
C.convertArrays2File([a,b], 'out.plt')
