# - usurp (array) -
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
m = C.addVar(m, 'ro')

c = C.array('cellN', ni-1, nj-1, 1)
c = C.initVars(c, 'cellN', 1)

import Post as P
r = P.usurp([m], [c]); print(r)

#=============================================================================
# Basic test case 2
#=============================================================================
m2 = G.cart((11,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.addVar(m2, 'ro')

c2 = C.array('cellN', ni-1, nj-1, 1)
c2 = C.initVars(c2, 'cellN', 1)
r = P.usurp([m,m2], [c,c2]); print(r)

#=============================================================================
# Basic test cas 3
#=============================================================================
m3 = G.cart((10,0,0), (10.*sqrt(2)/(ni-1),10.*sqrt(2)/(nj-1),1), (ni,nj,1))
m3 = T.rotate(m3, (10,0,0), (0,0,1), 45.)
m3 = C.addVar(m3, 'ro')

c3 = C.array('cellN', ni-1, nj-1, 1)
c3 = C.initVars(c3, 'cellN', 1)

C.convertArrays2File([m,m3], "out.plt", "bin_tp")

r = P.usurp([m,m3], [c,c3]); print(r)

m = T.addkplane(m)
m3 = T.addkplane(m3)

centers = C.node2Center([m,m3])
res1 = C.addVars([centers[0], r[0]])
res2 = C.addVars([centers[1], r[1]])
C.convertArrays2File([res1,res2], "centers.plt", "bin_tp")
