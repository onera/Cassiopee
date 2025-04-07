# - distribute (array) -
import Generator as G
import Distributor2 as D2
import KCore.test as test

# Premier cas pas trop gros
N = 11
arrays = []; off = 0
for i in range(N):
    a = G.cart( (off,0,0), (1,1,1), (10+i, 10, 10) )
    off += 9+i
    arrays.append(a)

out = D2.distribute(arrays, NProc=5, algorithm='gradient')
test.testO(out['distrib'], 1)

out = D2.distribute(arrays, NProc=5, algorithm='genetic')
test.testO(out['distrib'], 2)

out = D2.distribute(arrays, NProc=5, algorithm='fast')
test.testO(out['distrib'], 3)

# Deuxieme cas plus gros
N = 11
arrays = []; offj = 0
for j in range(N):
    offi = 0
    for i in range(N):
        a = G.cart( (offi,offj,0), (1,1,1), (10+i, 10+j, 10) )
        offi += 9+i
        arrays.append(a)
    offj += 9+j

out = D2.distribute(arrays, NProc=50, algorithm='gradient')
test.testO(out['distrib'], 4)

out = D2.distribute(arrays, NProc=50, algorithm='genetic')
test.testO(out['distrib'], 5)

out = D2.distribute(arrays, NProc=50, algorithm='fast')
test.testO(out['distrib'], 6)
