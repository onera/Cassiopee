# - computeDiv2 (pyTree) -
import Converter.PyTree as C
import Post.PyTree      as P
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test       as test

# --------
# 3D NGON
#--------
ni = 15; nj = 15; nk = 15
m1 = G.cartNGon(( 0, 0, 0), (10./(ni-1),10./(nj-1),2./(nk-1)), (ni,nj,nk))
m2 = G.cartNGon(( 0,10, 0), (10./(ni-1),10./(nj-1),2./(nk-1)), (ni,nj,nk))
m3 = G.cartNGon(( 0, 0, 2), (10./(ni-1),10./(nj-1),2./(nk-1)), (ni,nj,nk))

t  = C.newPyTree(['Base',m1,m2,m3])

t  = X.connectMatch(t,dim=3)
t  = C.fillEmptyBCWith(t,"wall",'BCWall')

C._initVars(t, '{centers:fldX}= ({centers:CoordinateX})**2')
C._initVars(t, '{centers:fldY}= ({centers:CoordinateY})**2')
C._initVars(t, '{centers:fldZ}= ({centers:CoordinateZ})**2')

C._initBCDataSet(t,'{fldX}=0.5')
C._initBCDataSet(t,'{fldY}=0.5')
C._initBCDataSet(t,'{fldZ}=0.5')

G._getVolumeMap(t)

P._computeDiv2(t, 'centers:fld')
test.testT(t,1)
