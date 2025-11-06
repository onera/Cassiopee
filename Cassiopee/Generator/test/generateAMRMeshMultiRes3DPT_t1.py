import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Geom.IBM as D_IBM
import Generator.AMR as G_AMR
import KCore.test as test

LOCAL      = test.getLocal()
snear      = 0.4
dFar       = 5
levelMax   = 2
vmins      = [[6,6]]
vminsTbox  = [[5,5]]
dim        = 3
dir_sym    = 2
isSymmetry = True

tb    = C.convertFile2PyTree('m6cleanDouble.cgns')
#tbox  = C.convertFile2PyTree('tboxDouble3D.cgns')
#D_IBM._setSnear(tbox,2*snear)
tbox = None
D_IBM._setSnear(tb,snear)
D_IBM._setDfar(tb,dFar)
D_IBM._symetrizePb(tb, 'Base', snear, dir_sym)
Cmpi.barrier()
opt   = True

t_AMR = G_AMR.generateAMRMesh(tb=tb, toffset=None, levelMax=levelMax, vmins=vmins, dim=dim,
                              check=False, opt=False, octreeMode=1, tbox=tbox, localDir=LOCAL,
                              NumMinDxLarge=0)
test.testT(t_AMR,1)
#Cmpi.convertPyTree2File(t_AMR,'check_multi_t1_3D.cgns')
