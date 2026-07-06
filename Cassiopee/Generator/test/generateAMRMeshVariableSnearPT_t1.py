# - generateAMRMesh (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.AMR as G_AMR
import Geom.PyTree as D
import Geom.IBM as D_IBM
import Transform.PyTree as T
import KCore.test as test

LOCAL = test.getLocal()

# 2D
a = D.naca(12.)
dimPb = 2
D_IBM._setSnear(a, 0.002)
D_IBM._setDfar(a, 10.)

a = T.splitNParts(a, 3, dirs=[1])
for i, z in enumerate(Internal.getZones(a)):
    D_IBM._setSnear(z, 0.002*(i+1))
tb = C.newPyTree(["BODY",a])
t = G_AMR.generateAMRMesh(tb, vmins=[[5,5]], dim=dimPb, octreeMode=1, check=False, localDir=LOCAL)
test.testT(t,1)
#C.convertPyTree2File(t,'check_t1_2D.cgns')
