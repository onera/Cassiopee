# - printStats (pyTree) -
import Generator.PyTree as G
import Distributor2.PyTree as D2
import Converter.PyTree as C
import Connector.PyTree as X

t = C.newPyTree(['Base'])
pos = 0
for i in range(11):
    a = G.cart((pos,0,0), (1,1,1), (10+i, 10, 10))
    pos += 10 + i - 1
    t[2][1][2].append(a)
t = X.connectMatch(t)

# Distribute on 3 processors
t, stats = D2.distribute(t, 3)

D2.printStats(t)
#>> Info: varMin=3.636364%, varMax=12.727273%, varRMS=5.352582%
#>> Info: external com ratio=20.000000%
