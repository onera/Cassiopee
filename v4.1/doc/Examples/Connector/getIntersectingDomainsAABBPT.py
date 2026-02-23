# - getIntersectingDomainsAABB (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Connector.PyTree as X
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((9.,0,0), (1,1,1), (10,10,10))

Ni = 50; Nj = 50; Nk = 2
a = G.cart((0,0,0),(1./(Ni-1), 1./(Nj-1),1), (Ni,Nj,Nk))
b = G.cart((0,0,0),(2./(Ni-1), 2./(Nj-1),1), (Ni,Nj,Nk)); b[0] = 'cart2'
a = T.rotate(a, (0,0,0), (0,0,1), 10.); a = T.translate(a, (0.5,0.5,0))

bb = G.BB([a,b])
ret = X.getIntersectingDomainsAABB(bb)
print(ret)
