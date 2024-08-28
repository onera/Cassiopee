# - cart (array) -
import Converter as C
import Generator as G
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
C.convertArrays2File(a, 'out.plt')
