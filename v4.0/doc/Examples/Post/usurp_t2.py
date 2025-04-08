# - usurp -
import Post as P
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test
from math import *

#=============================================================================
# Basic test case 1
#=============================================================================
ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m1 = C.addVars(m1, 'ro')
c1 = C.array('cellN', ni-1, nj-1, 1)
c1 = C.initVars(c1, 'cellN', 1)

r = P.usurp([m1], [c1])
if r is not None: test.testA(r,1)

#=============================================================================
# Basic test case 2
#=============================================================================
m2 = G.cart((11,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.addVars(m2, 'ro')
c2 = C.array('cellN', ni-1, nj-1, 1)
c2 = C.initVars(c2, 'cellN', 1)
r = P.usurp([m1,m2], [c1,c2])
if r is not None: test.testA(r,2)

#=============================================================================
# Basic test cas 3
#=============================================================================
m3 = G.cart((10,0,0), (10.*sqrt(2)/(ni-1),10.*sqrt(2)/(nj-1),1), (ni,nj,1))
m3 = T.rotate(m3, (10,0,0), (0,0,1), 45.)
m3 = C.addVars(m3, 'ro')

c3 = C.array('cellN', ni-1, nj-1, 1)
c3 = C.initVars(c3, 'cellN', 1)
r = P.usurp([m1,m3], [c1,c3])
if r is not None: test.testA([r[1]],3)
