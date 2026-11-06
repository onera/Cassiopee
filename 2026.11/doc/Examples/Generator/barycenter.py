# - barycenter (array) -
import Generator as G
import Converter as C
a = G.cart((0.,0.,0.), (0.1,0.1,1.), (20,20,20))
print(G.barycenter(a))
#>> [0.9500000000000001, 0.9500000000000005, 9.5]
w = C.initVars(a, 'weight', 1.); w = C.extractVars(w,['weight'])
print(G.barycenter(a, w))
#>> [0.9500000000000001, 0.9500000000000005, 9.5]