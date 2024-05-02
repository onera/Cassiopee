# - getNobNozOfZone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base', 'Base2'])
t[2][1][2] += [a]; t[2][2][2] += [b]
(nob, noz) = C.getNobNozOfZone(a, t); print(nob, noz)
#>> 1 0
# This means that t[2][nob][2][noz] = a
