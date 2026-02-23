# - computeDiv2 (pyTree) -
import Converter.PyTree as C
import Post.PyTree      as P
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test       as test

def F(x): return 2*x

#---------------
# 2D STRUCT - XY
#---------------
ni = 15; nj = 15
m1 = G.cart(( 0, 0, 0), (10./(ni-1),10./(nj-1),1.), (ni,nj,1))
m2 = G.cart(( 0,10, 0), (10./(ni-1),10./(nj-1),1.), (ni,nj,1))
m3 = G.cart((10, 0, 0), (10./(ni-1),10./(nj-1),1.), (ni,nj,1))

t  = C.newPyTree(['Base',m1,m2,m3])

t  = X.connectMatch(t,dim=2)
t  = C.fillEmptyBCWith(t,"wall",'BCWall')

C._initVars(t, '{centers:fld1X}= ({centers:CoordinateX})**2')
C._initVars(t, '{centers:fld1Y}= ({centers:CoordinateY})**2')

C._initBCDataSet(t,'{fld1X}=0.5')
C._initBCDataSet(t,'{fld1Y}=0.5')

C._initVars(t, 'centers:fld2X', F, ['centers:fld1X'])
C._initVars(t, 'centers:fld2Y', F, ['centers:fld1Y'])

C._initBCDataSet(t,'{fld2X}=1.')
C._initBCDataSet(t,'{fld2Y}=1.')

P._computeDiv2(t, ['centers:fld1', 'centers:fld2'], rmVar=True)
test.testT(t,1)


#---------------
# 2D STRUCT - YZ
#---------------
m1 = G.cart(( 0, 0, 0), (1.,10./(ni-1),10./(nj-1)), (1,ni,nj))
m2 = G.cart(( 0, 0,10), (1.,10./(ni-1),10./(nj-1)), (1,ni,nj))
m3 = G.cart(( 0,10, 0), (1.,10./(ni-1),10./(nj-1)), (1,ni,nj))

t  = C.newPyTree(['Base',m1,m2,m3])

t  = X.connectMatch(t,dim=2)
t  = C.fillEmptyBCWith(t,"wall",'BCWall')

C._initVars(t, '{centers:fld1Y}= ({centers:CoordinateY})**2')
C._initVars(t, '{centers:fld1Z}= ({centers:CoordinateZ})**2')

C._initBCDataSet(t,'{fld1Y}=0.5')
C._initBCDataSet(t,'{fld1Z}=0.5')

C._initVars(t, 'centers:fld2Y', F, ['centers:fld1Y'])
C._initVars(t, 'centers:fld2Z', F, ['centers:fld1Z'])

C._initBCDataSet(t,'{fld2Y}=1.')
C._initBCDataSet(t,'{fld2Z}=1.')

P._computeDiv2(t, ['centers:fld1', 'centers:fld2'], rmVar=True)
test.testT(t,2)


#---------------
# 2D STRUCT - XZ
#---------------
m1 = G.cart(( 0, 0, 0), (10./(ni-1),1.,10./(nj-1)), (ni,1,nj))
m2 = G.cart(( 0, 0,10), (10./(ni-1),1.,10./(nj-1)), (ni,1,nj))
m3 = G.cart((10, 0, 0), (10./(ni-1),1.,10./(nj-1)), (ni,1,nj))

t  = C.newPyTree(['Base',m1,m2,m3])

t  = X.connectMatch(t,dim=2)
t  = C.fillEmptyBCWith(t,"wall",'BCWall')

C._initVars(t, '{centers:fld1X}= ({centers:CoordinateX})**2')
C._initVars(t, '{centers:fld1Z}= ({centers:CoordinateZ})**2')

C._initBCDataSet(t,'{fld1X}=0.5')
C._initBCDataSet(t,'{fld1Z}=0.5')

C._initVars(t, 'centers:fld2X', F, ['centers:fld1X'])
C._initVars(t, 'centers:fld2Z', F, ['centers:fld1Z'])

C._initBCDataSet(t,'{fld2X}=1.')
C._initBCDataSet(t,'{fld2Z}=1.')

P._computeDiv2(t, ['centers:fld1', 'centers:fld2'], rmVar=True)
test.testT(t,3)


#----------
# 3D STRUCT
#----------
ni = 15; nj = 15; nk = 15
m1 = G.cart(( 0, 0, 0), (10./(ni-1),10./(nj-1),2./(nk-1)), (ni,nj,nk))
m2 = G.cart(( 0,10, 0), (10./(ni-1),10./(nj-1),2./(nk-1)), (ni,nj,nk))
m3 = G.cart(( 0, 0, 2), (10./(ni-1),10./(nj-1),2./(nk-1)), (ni,nj,nk))

t  = C.newPyTree(['Base',m1,m2,m3])
t  = X.connectMatch(t,dim=3)
t  = C.fillEmptyBCWith(t,"wall",'BCWall')

C._initVars(t, '{centers:fld1X}= ({centers:CoordinateX})**2')
C._initVars(t, '{centers:fld1Y}= ({centers:CoordinateY})**2')
C._initVars(t, '{centers:fld1Z}= ({centers:CoordinateZ})**2')

C._initBCDataSet(t,'{fld1X}=0.5')
C._initBCDataSet(t,'{fld1Y}=0.5')
C._initBCDataSet(t,'{fld1Z}=0.5')

C._initVars(t, 'centers:fld2X', F, ['centers:fld1X'])
C._initVars(t, 'centers:fld2Y', F, ['centers:fld1Y'])
C._initVars(t, 'centers:fld2Z', F, ['centers:fld1Z'])

C._initBCDataSet(t,'{fld2X}=1.')
C._initBCDataSet(t,'{fld2Y}=1.')
C._initBCDataSet(t,'{fld2Z}=1.')

P._computeDiv2(t, ['centers:fld1', 'centers:fld2'], rmVar=True)
test.testT(t,4)


#--------
# 3D NGON
#--------
ni = 30; nj = 40; nk = 10
t = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(t, '{centers:fld1X}= cos({centers:CoordinateX})')
C._initVars(t, '{centers:fld1Y}= 4.*{centers:CoordinateY}')
C._initVars(t, '{centers:fld1Z}= {centers:CoordinateY}*{centers:CoordinateZ}**2.')

C._initVars(t, 'centers:fld2X', F, ['centers:fld1X'])
C._initVars(t, 'centers:fld2Y', F, ['centers:fld1Y'])
C._initVars(t, 'centers:fld2Z', F, ['centers:fld1Z'])

P._computeDiv2(t, ['centers:fld1', 'centers:fld2'], rmVar=True)
test.testT(t,5)
