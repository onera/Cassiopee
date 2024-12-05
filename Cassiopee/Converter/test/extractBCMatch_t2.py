# - extractBCMatch (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as CP
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

a = G.cart((1,1,1), (1.,1.,1.), (4,10,3)); a[0]='cart1'
b = G.cart((4,2,0), (1.,1.,1.), (5,8,5)) ; b[0]='cart2'
c = G.cart((4,9,1), (1.,1.,1.), (4,5,6)); c[0]='cart3'

a = T.reorder(a,(1,2,3))
b = T.reorder(b,(2,-1,3))
c = T.reorder(c,(3,2,-1))

t = CP.newPyTree(['Base',a,b,c])

t = CP.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = CP.initVars(t, '{centers:G}=2.3')
t = CP.initVars(t, '{centers:H}={centers:CoordinateY}')
t = CP.initVars(t, '{centers:M}={centers:CoordinateX}')
t = X.connectMatch(t,dim=3)
t = CP.fillEmptyBCWith(t,"wall",'BCWall')

zones = Internal.getZones(t)
it    = 0

for z in zones:
    dim = Internal.getZoneDim(z)
    gcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')

    for gc in gcs:
        zname  = Internal.getValue(gc)
        zdonor = Internal.getNodeFromName(t,zname)

        [indFaceR,fldFace]  = CP.extractBCMatch(zdonor,gc,dim,['centers:G','centers:H','centers:M'])

        test.testO([indFaceR,fldFace], it)
        it  = it + 1

        [indFaceR,fldFace] = CP.extractBCMatch(zdonor,gc,dim)

        test.testO([indFaceR,fldFace], it)
        it  = it + 1
