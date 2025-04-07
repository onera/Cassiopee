import Converter as C
import Generator as G
import Transform as T
from math import *

#=============================================================================
# Basic test case 1
#=============================================================================
ni = 30; nj = 40
# premier maillage
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m1 = C.addVar(m, 'ro'); del m

c = C.array('cellN', ni-1, nj-1, 1)
c1 = C.initVars(c, 'cellN', 1); del c
#
import Post as P
r = P.usurp([m1], [c1])
print(r)

res = P.integ([m1],[c1],[])

#=============================================================================
# Basic test case 2
#=============================================================================
m = G.cart((11,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.addVar(m, 'ro'); del m

c = C.array('cellN', ni-1, nj-1, 1)
c2 = C.initVars(c, 'cellN', 1); del c
r = P.usurp([m1,m2], [c1,c2])
print(r)

#=============================================================================
# Basic test case 3
#=============================================================================
m = G.cart((10,0,0), (10.*sqrt(2)/(ni-1),10.*sqrt(2)/(nj-1),1), (ni,nj,1))
mp = T.rotate(m, (10,0,0), (0,0,1), 45.); del m
m3 = C.addVar(mp, 'ro'); del mp

c = C.array('cellN', ni-1, nj-1, 1)
c3 = C.initVars(c, 'cellN', 1); del c

C.convertArrays2File([m1,m3], "out.plt", "bin_tp")

r = P.usurp([m1,m3], [c1,c3])
print(r)

m4 = T.addkplane(m1); del m1
m5 = T.addkplane(m3); del m3

c1 = P.node2Center(m4)
c2 = P.node2Center(m5)
res1 = C.addVars([c1, r[0]])
res2 = C.addVars([c2, r[1]])
C.convertArrays2File([res1,res2], "centers.plt", "bin_tp")
