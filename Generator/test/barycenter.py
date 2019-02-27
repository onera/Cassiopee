# - barycenter (array) -
import Generator as G
import Converter as C
a = G.cart((0.,0.,0.), (0.1,0.1,1.), (20,20,20))
print(G.barycenter(a))
w = C.initVars(a, 'weight', 1.); w = C.extractVars(w,['weight'])
print(G.barycenter(a, w))
