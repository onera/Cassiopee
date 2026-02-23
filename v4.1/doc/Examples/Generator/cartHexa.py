# - cartHexa (array) -
import Generator as G
import Converter as C

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
C.convertArrays2File([a], 'out.plt')
