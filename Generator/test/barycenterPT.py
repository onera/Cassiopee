# - barycenter (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
a = G.cart((0.,0.,0.), (0.1,0.1,1.), (20,20,20))
print G.barycenter(a)
a = C.initVars(a, 'weight', 1.)
print G.barycenter(a, 'weight')
