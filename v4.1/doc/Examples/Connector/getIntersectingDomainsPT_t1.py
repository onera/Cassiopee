# - getIntersectingDomains (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import KCore.test as test

t = C.newPyTree(['Base1'])
Ni = 4; Nj = 4; Nk = 4; dx = 0.
for i in range(4):
    z = G.cart((dx,dx,dx),(1./(Ni-1), 1./(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
    t[2][1][2] += [z]
    dx += 0.3

interDict = X.getIntersectingDomains(t, method='AABB')
test.testO(interDict,1)

interDict = X.getIntersectingDomains(t, method='OBB')
test.testO(interDict,2)

interDict = X.getIntersectingDomains(t, method='hybrid')
test.testO(interDict,3)
