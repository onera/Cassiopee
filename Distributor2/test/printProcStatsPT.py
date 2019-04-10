# - printProcStats (pyTree) -
import Generator.PyTree as G
import Distributor2.PyTree as D2
import Converter.PyTree as C
import Connector.PyTree as X

N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in range(N):
    a = G.cart((pos,0,0), (1,1,1), (10+i, 10, 10))
    pos += 10 + i - 1
    t[2][1][2].append(a)
t = X.connectMatch(t)

# Distribute on 3 processors
t, stats = D2.distribute(t, 3)

# With stats and NProc
D2.printProcStats(t, stats, NProc=3)
# NProc is guessed from stats
D2.printProcStats(t, stats)
# All are guessed from t
D2.printProcStats(t)

C.convertPyTree2File(t, 'out.cgns')
