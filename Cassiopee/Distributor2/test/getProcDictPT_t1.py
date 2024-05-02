# - getProcDict (pyTree) -
import Generator.PyTree as G
import Distributor2.PyTree as D2
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in range(N):
    a = G.cart( (pos,0,0), (1,1,1), (10+i, 10, 10) )
    pos += 10 + i - 1
    t[2][1][2].append(a)

t = X.connectMatch(t)
t, stats = D2.distribute(t, 3)

proc = D2.getProcDict(t)
test.testO(proc, 1)

proc = D2.getProcDict(t, prefixByBase=True)
test.testO(proc, 2)
