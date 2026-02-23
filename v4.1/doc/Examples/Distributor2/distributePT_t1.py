# - distribute (pyTree) -
import Generator.PyTree as G
import Distributor2.PyTree as D2
import Converter.PyTree as C
import Connector.PyTree as X
import KCore.test as test

N = 11

# Test raccords matchs
t = C.newPyTree(['Base'])
off = 0
for i in range(N):
    a = G.cart( (off,0,0), (1,1,1), (10+i, 10, 10) )
    off += 9+i
    t[2][1][2].append(a)
t = X.connectMatch(t)

t, stats = D2.distribute(t, NProc=5, algorithm='gradient')
test.testT(t, 1)

t, stats = D2.distribute(t, NProc=5, algorithm='genetic')
test.testT(t, 2)

t, stats = D2.distribute(t, NProc=5, algorithm='fast')
test.testT(t, 3)

# avec overlap
t = C.newPyTree(['Base'])
off = 0
for i in range(N):
    a = G.cart( (off,0,0), (1,1,1), (10+i, 10, 10) )
    C._addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
    C._addBC2Zone(a, 'overlap', 'BCOverlap', 'imax')
    off += 9+i
    t[2][1][2].append(a)

t, stats = D2.distribute(t, NProc=5, algorithm='gradient')
test.testT(t, 4)

t, stats = D2.distribute(t, NProc=5, algorithm='genetic')
test.testT(t, 5)

t, stats = D2.distribute(t, NProc=5, algorithm='fast')
test.testT(t, 6)

# Avec zones qui se recouvrent
t = C.newPyTree(['Base'])
off = 0
for i in range(N):
    a = G.cart( (off,0,0), (1,1,1), (10+i, 10, 10) )
    off += 9+i
    t[2][1][2].append(a)

t, stats = D2.distribute(t, NProc=5, algorithm='gradient', useCom='bbox')
test.testT(t, 7)

# Avec des poids
t = C.newPyTree(['Base'])
off = 0
weightDict={}
for i in range(N):
    a = G.cart( (off,0,0), (1,1,1), (10+i, 10, 10) )
    off += 9+i
    t[2][1][2].append(a)
    weightDict[a[0]]=i+1

t, stats = D2.distribute(t, NProc=5, weight=weightDict,algorithm='gradient', useCom='bbox')
test.testT(t,8)

t2 = C.convertArray2Hexa(t)
stats = D2._distribute(t2, NProc=5, weight=weightDict,algorithm='gradient', useCom='bbox', mode='cells')
test.testT(t2,9)

t2 = C.convertArray2NGon(t)
stats = D2._distribute(t2, NProc=5, weight=weightDict,algorithm='gradient', useCom='bbox', mode='cells')
test.testT(t2,10)
