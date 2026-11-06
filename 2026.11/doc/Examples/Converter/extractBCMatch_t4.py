# - extractBCMatch (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as CP
import Converter.Internal as Internal
import Connector.PyTree as X
import KCore.test as test

a = G.cartNGon((1,1,1), (1.,1.,1.), (4,10,1)); a[0]='cart1'
b = G.cartNGon((4,2,1), (1.,1.,1.), (5,8,1)) ; b[0]='cart2'
c = G.cartNGon((4,9,1), (1.,1.,1.), (4,5,1)); c[0]='cart3'

t = CP.newPyTree(['Base',a,b,c])

t = CP.initVars(t, '{F}=3*{CoordinateX}+2*{CoordinateY}')
t = CP.initVars(t, '{centers:G}=2.3')
t = CP.initVars(t, '{centers:H}={centers:CoordinateY}')
t = CP.initVars(t, '{centers:M}={centers:CoordinateX}')
t = X.connectMatch(t,dim=2)
t = CP.fillEmptyBCWith(t,"wall",'BCWall')

varL = ['H']

zones = Internal.getZones(t)
it    = 0

for z in zones:
    dim = Internal.getZoneDim(z)
    gcs = Internal.getNodesFromType2(z, 'GridConnectivity_t')

    for gc in gcs:
        it  = it + 1
        zname  = Internal.getValue(gc)
        zdonor = Internal.getNodeFromName(t,zname)

        [indFaceR,fldFace]  = CP.extractBCMatch(zdonor,gc,dim,varL)
        test.testO([indFaceR,fldFace], it)
