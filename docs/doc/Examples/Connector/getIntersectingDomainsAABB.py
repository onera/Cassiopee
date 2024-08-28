# - getIntersectingDomainsAABB (array) -
import Generator as G
import Converter as C
import Connector as X

a = G.cart((0.,0,0), (1,1,1), (10,10,10))
b = G.cart((9.,0,0), (1,1,1), (10,10,10))

bb = G.BB([a,b])
ret = X.getIntersectingDomainsAABB(bb)
print(ret)
