# - setNANValuesAt (array) -
import Generator as G
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 1.)
a = C.setNANValuesAt(a, 'F', 0.)
