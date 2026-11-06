# - integ (array) -
import Converter as C
import Generator as G
import Transform as T
import Post as P

from math import *

#=============================================================================
# Integration coord noeuds, F centres, surfaces planes
#=============================================================================
ni = 30; nj = 40
# premier maillage (k=1)
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.addVars(m, 'ro')
c = C.array('vx,vy,vz', ni-1, nj-1, 1)
c = C.initVars(c, 'vx,vy,vz', 1.)
res1 = P.integ([m], [c], [])
print('integ scal:', res1)
resn1 = P.integNorm([m], [c], [])
print('integ norm:', resn1)

# Deuxieme maillage (j=1)
m = G.cart((0,0,0), (10./(ni-1),1,10./(nj-1)), (ni,1,nj))
c = C.array('vx,vy,vz', ni-1, 1, nj-1)
c = C.initVars(c, 'vx,vy,vz', 1.)
res2 = P.integ([m], [c], [])
print('integ scal:', res2)
resn2 = P.integNorm([m], [c], [])
print('integ norm:', resn2)

# Troisieme maillage (i=1)
m3 = G.cart((0,0,0), (1, 10./(ni-1),10./(nj-1)), (1,ni,nj))
c = C.array('vx,vy,vz', 1, ni-1, nj-1)
c3 = C.initVars(c, 'vx,vy,vz', 1.); del c
res3 = P.integ([m3], [c3], [])
print('integ scal:', res3)
resn3 = P.integNorm([m3], [c3], [])
print('integ norm:', resn3)

#=============================================================================
# Integration coord noeuds, F centres, surfaces courbes
#=============================================================================
ni = 30; nk = 40
m = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (ni,2,nk))
m1 = T.subzone(m, (1,2,1), (ni,2,nk))
#C.convertArrays2File( [m1], "out.plt", "bin_tp")
c = C.array('vx,vy,vz', ni-1, 1, nk-1)
c1 = C.initVars(c, 'vx,vy,vz', 1.); del c
res1 = P.integ([m1], [c1], [])
print('integ scal:', res1, pi*2*1.*10)
resn1 = P.integNorm([m1], [c1], [])
print('integ norm:', resn1)

#res3 = P.integNormProduct([m1],[c1],[])
#print(res3)
#print(res1, res2, res3)
