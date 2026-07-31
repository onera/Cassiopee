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
snear = 0.002
D_IBM._setSnear(a, snear)
D_IBM._setDfar(a, 10.)

a = T.splitNParts(a, 3, dirs=[1])
for i, z in enumerate(Internal.getZones(a)):
    D_IBM._setSnear(z, snear*(i+1))
tb = C.newPyTree(["BODY",a])


# Tboxs
l1 = D.line((1.55,-0.20, 0.0), (2.30,-0.20, 0.0), 200)
l2 = D.line((2.30,-0.20, 0.0), (2.30, 0.85, 0.0), 200)
l3 = D.line((2.30, 0.85, 0.0), (1.55, 0.85, 0.0), 200)
l4 = D.line((1.55, 0.85, 0.0), (1.55,-0.20, 0.0), 200)
a  = T.join([l1,l2,l3,l4])
D_IBM._setSnear(a, snear*8)

b  = T.translate(a, (0.,-1.30,0.))
D_IBM._setSnear(b, snear*16)
tbox = C.newPyTree(["BODY1",a, 'COPY1', b])
D_IBM._setDfar(tbox,10)
D_IBM._setIBCType(tbox,"None")

t = G_AMR.generateAMRMesh(tb, vmins=[[5,5]], dim=dimPb, octreeMode=1, tbox=tbox, check=False, localDir=LOCAL)
test.testT(t,1)
#C.convertPyTree2File(t,'check_t1_2D.cgns')
