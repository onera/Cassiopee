# - cartNGon (array) -
import Generator as G
import Converter as C

a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.2), (20,20,20))
C.convertArrays2File(a, 'out.plt')
