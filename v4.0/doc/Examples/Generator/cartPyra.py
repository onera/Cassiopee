# - cartHexa (array) -
import Generator as G
import Converter as C

a = G.cartPyra((0.,0.,0.), (1,1,1), (20,20,20))
C.convertArrays2File(a, 'out.tp')
