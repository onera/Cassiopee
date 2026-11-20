# - generateAMRMesh (pyTree) -
import Converter.PyTree as C
import Generator.AMR as G_AMR
import Generator.PyTree as G
import Geom.PyTree as D
import Post.PyTree as P
import Transform.PyTree as T
import Converter.Internal as Internal
import Geom.IBM as D_IBM
import KCore.test as test

LOCAL = test.getLocal()

tb1 = G.cart((0.,0.,0.), (0.1,0.1,0.1), (10,20,5))
tb1 = P.exteriorFaces(tb1)
tb1 = T.splitSharpEdges(tb1, 30.)
zonesJoin = []
for z in Internal.getZones(tb1):
    maxval = C.getMaxValue(z, 'CoordinateY');
    if maxval > 1e-06:
        zonesJoin.append(z)
        Internal._rmNode(tb1, z)
    else:
        z[0]="midPlane"
zonesJoin = T.join(zonesJoin)
tb = C.newPyTree(['Base', tb1, zonesJoin])

dictExtension = {
    "midPlane": -1
}
dictTolerance= {
    "midPlane": 1e-05
}
G_AMR._addExtensionInfo(tb, dictExtension, dictTolerance=dictTolerance)
test.testT(tb,1)

tbExtSym = G_AMR.createExtension__(tb)
test.testT(tbExtSym,2)
