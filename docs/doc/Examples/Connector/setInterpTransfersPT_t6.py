import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test

d = 2
a = D.sphere6((0,0,0),0.5,N=20)
dk = G.cart((0,0,0),(0.01,1,1),(11,1,1))
a = G.addNormalLayers(a,dk)
t = C.newPyTree(['Base']); t[2][1][2] = a
t = X.connectMatch(t,dim=3)
Internal._addGhostCells(t,t,d,adaptBCs=1)
#---------
# Centers
#---------
C._initVars(t,'centers:F',0.)
tc = C.node2Center(t)
C._initVars(tc,'{F}={CoordinateX}*{CoordinateY}')
# stockage direct
t1 = X.setInterpData(t, tc, storage='direct', loc='centers',itype='abutting')
X._setInterpTransfers(t1,tc,variables=['F'])
test.testT(t1,1)

# stockage inverse
tc1 = X.setInterpData(t, tc, storage='inverse', loc='centers',itype='abutting')
t1 = X.setInterpTransfers(t,tc1,variables=['F'])
test.testT(t1,2)
#---------
# Nodes
#---------
# noeuds, stockage direct
C._rmVars(t,['centers:F'])
C._initVars(t,'F', 0.)
td = C.initVars(t,'{F}={CoordinateX}*{CoordinateY}')
t1 = X.setInterpData(t, td, storage='direct', loc='nodes',itype='abutting')
X._setInterpTransfers(t1,td,variables=['F'])
test.testT(t1,3)

# noeuds, stockage inverse
td1 = X.setInterpData(t, td, storage='inverse', loc='nodes',itype='abutting')
t1 = X.setInterpTransfers(t,td1,variables=['F'])
test.testT(t1,4)
