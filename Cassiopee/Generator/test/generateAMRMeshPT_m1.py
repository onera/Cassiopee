
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Generator.PyTree as G
import Geom.PyTree as D
import Geom.IBM as D_IBM
import KCore.test as test

#2D
a = D.naca(12.)
dimPb = 2
D_IBM._setSnear(a, 0.5)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 10.)

tb = C.newPyTree(["BODY",a])
Cmpi.barrier()

t = G.generateAMRMesh(tb, levelMax=3, vmins=[15,15], dim=dimPb, check=False)
if Cmpi.rank == 0: test.testT(t,1)
else: test.testT(t,12)

# 3D
a = D.sphere((0.,0.,0.),0.1)
dimPb = 3
D_IBM._setSnear(a, 0.1)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 5.)

tb = C.newPyTree(["BODY",a])
Cmpi.barrier()


tb = C.newPyTree(["BODY",a])
toffset = C.newPyTree(['R1'])
toffset[2][1][2] = [D.sphere((0.,0.,0.),0.5)]
t = G.generateAMRMesh(tb, toffset=toffset, levelMax=4, vmins=[5], dim=dimPb, check=False)
if Cmpi.rank == 0: test.testT(t,2)
else: test.testT(t,22)
