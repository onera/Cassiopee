
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Generator.PyTree as G
import Geom.PyTree as D
import Geom.IBM as D_IBM
import Connector.AMR as X_AMR

a = D.naca(12.)
dimPb = 2
D_IBM._setSnear(a, 0.05)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 20.)

tb = C.newPyTree(["BODY",a])
if Cmpi.rank==0: C.convertPyTree2File(tb,"case.cgns")
Cmpi.barrier()
t = G.generateAMRMesh(tb, levelMax=6, vmins=[15,10,10,8,5], dim=dimPb, check=True)
Cmpi.convertPyTree2File(t,"out.cgns")
