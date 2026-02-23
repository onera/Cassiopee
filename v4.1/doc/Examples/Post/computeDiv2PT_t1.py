# - computeDiv2 (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
C._initVars(m, '{centers:fldX}= cos({centers:CoordinateX})')
C._initVars(m, '{centers:fldY}= 4.*{centers:CoordinateY}')
C._initVars(m, '{centers:fldZ}= {centers:CoordinateY}*{centers:CoordinateZ}**2.')
P._computeDiv2(m, 'centers:fld')
test.testT(m,1)

ni = 30; nj = 40; nk = 10
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(m, '{centers:fldX}= cos({centers:CoordinateX})')
C._initVars(m, '{centers:fldY}= 4.*{centers:CoordinateY}')
C._initVars(m, '{centers:fldZ}= {centers:CoordinateY}*{centers:CoordinateZ}**2.')
P._computeDiv2(m, 'centers:fld')
test.testT(m,2)

#-----
# NGON
#-----
ni = 30; nj = 40
m = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
C._initVars(m, '{centers:fldX}= cos({centers:CoordinateX})')
C._initVars(m, '{centers:fldY}= 4.*{centers:CoordinateY}')
C._initVars(m, '{centers:fldZ}= {centers:CoordinateY}*{centers:CoordinateZ}**2.')
P._computeDiv2(m, 'centers:fld')
test.testT(m,3)

ni = 30; nj = 40; nk = 10
m = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(m, '{centers:fldX}= cos({centers:CoordinateX})')
C._initVars(m, '{centers:fldY}= 4.*{centers:CoordinateY}')
C._initVars(m, '{centers:fldZ}= {centers:CoordinateY}*{centers:CoordinateZ}**2.')
P._computeDiv2(m, 'centers:fld')
test.testT(m,4)
