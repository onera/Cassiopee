# join avec les coordonnees et champs desordonnes
import Converter.PyTree as C
import Transform.PyTree as T
import Geom.PyTree as D
import Converter.Internal as Internal
import KCore.test as test

a = D.sphere6((0.,0.,0.),1.,N=10)
C._initVars(a,'Fx=1.')
C._initVars(a,'Fy=2.')
C._initVars(a,'Fz=3.')
C._initVars(a[1],'toto=0.')

noz = 0
for z in Internal.getZones(a):
    if noz != 1:
        C._initVars(z,'{centers:Gx}=2.')
        C._initVars(z,'{centers:Gy}=1.')
    else:
        C._initVars(z,'{centers:Gy}=1.')
        C._initVars(z,'{centers:Gx}=2.')
    GC = Internal.getNodeFromType(z,'GridCoordinates_t')
    Internal._rmNodesFromName(z,GC[0])
    if noz == 1:
        XA = Internal.getNodeFromName(GC,'CoordinateX')
        Internal._rmNodesFromName(GC,'CoordinateX')
        GC[2].append(XA)
    z[2].append(GC)
    noz+=1
res = T.join(a[0:2])
test.testT(res,1)
