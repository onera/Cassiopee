# - splitNParts (array) -
import Generator as G
import Transform as T
import Converter as C
import KCore.test as test
a = G.cart((0,0,0), (1,1,1), (41,41,41))
b = G.cart((40,0,0), (1,1,1), (21,41,21))
[a,b] = C.initVars([a,b], 'f', 1.)
N = 10
c = 1
for mg in range(3):
    for dirs0 in [[1],[2],[3],[1,2],[1,3],[2,3],[1,2,3]]:
        res = T.splitNParts([a,b], N, multigrid=mg, dirs=dirs0)
        test.testA(res, c)
        c += 1
