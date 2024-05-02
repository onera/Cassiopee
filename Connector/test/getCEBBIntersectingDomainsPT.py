# - getCEBBIntersectingDomains (pyTree) -
import Connector.PyTree as X
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10)); a[0] = 'cart1'
b = G.cart((0.5,0.,0.),(0.1,0.1,0.1),(10,10,10)); b[0] = 'cart2'
c = G.cart((0.75,0.,0.),(0.1,0.1,0.1),(10,10,10)); c[0] = 'cart3'

t = C.newPyTree(['Cart']); t[2][1][2] += [a, b, c]
bases = Internal.getNodesFromType(t, 'CGNSBase_t')
base = bases[0]
doms = X.getCEBBIntersectingDomains(base, bases, 1); print(doms)
