# - getIntersectingDomains (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C

t = C.newPyTree(['Base1'])
Ni = 4; Nj = 4; Nk = 4; dx = 0.
for i in range(10):
    z = G.cart((dx,dx,dx),(1./(Ni-1), 1./(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
    t[2][1][2] += [z]
    dx += 0.3

interDict = X.getIntersectingDomains(t, method='hybrid')
print('Does cart.1 intersect cart.2 ?','cart.1' in interDict['cart.2'])
print('List of zones intersecting cart.2:', interDict['cart.2'])
