# - generateAMRMesh (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Generator.AMR as G_AMR
import Geom.PyTree as D
import Geom.IBM as D_IBM
import KCore.test as test

LOCAL = test.getLocal()

# 2D
a = D.naca(12.)
dimPb = 2
D_IBM._setSnear(a, 0.125)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 20.)

tb = C.newPyTree(["BODY",a])
t = G_AMR.generateAMRMesh(tb, levelMax=3, vmins=[[5,5]], dim=dimPb, check=False, localDir=LOCAL)
test.testT(t,1)
#C.convertPyTree2File(t,'check_t1_2D.cgns')

# 3D
a = D.sphere((0.,0.,0.),0.1)
dimPb = 3
D_IBM._setSnear(a, 0.5)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 5.)

tb = C.newPyTree(["BODY",a])
toffset = C.newPyTree(['R1'])
toffset[2][1][2] = [D.sphere((0.,0.,0.),0.5)]
D_IBM._setSnear(toffset, 0.5)
t = G_AMR.generateAMRMesh(tb, toffset=toffset, levelMax=2, vmins=[[5]], dim=dimPb, check=False, localDir=LOCAL)
test.testT(t,2)
#C.convertPyTree2File(t,'check_t1_3D.cgns')
#C.convertPyTree2File(toffset,'check_t1_3D_offset.cgns')
#C.convertPyTree2File(a,'check_t1_3D_sphere.cgns')
