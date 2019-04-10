# - copyDistribution (pyTree) -
import Converter.PyTree as C
import Distributor2.PyTree as D2
import Generator.PyTree as G

# Case
N = 11
t = C.newPyTree(['Base'])
pos = 0
for i in range(N):
    a = G.cart((pos,0,0), (1,1,1), (10+i, 10, 10))
    a[0] = 'cart%d'%i
    pos += 10 + i - 1
    D2._addProcNode(a, i)
    t[2][1][2].append(a)

t2 = C.newPyTree(['Base'])
for i in range(N):
    a = G.cart((pos,0,0), (1,1,1), (10+i, 10, 10))
    a[0] = 'cart%d'%i
    pos += 10 + i - 1
    t2[2][1][2].append(a)
t2 = D2.copyDistribution(t2, t)
C.convertPyTree2File(t2, 'out.cgns')
