import Converter.PyTree as C
import Geom.IBM as D_IBM
import Generator.AMR as G_AMR
import KCore.test as test

LOCAL      = test.getLocal()
snear      = 0.001
dFar       = 10
levelMax   = 20
vmins      = [15,14,14,10,11,5,4]
vminsTbox  = [[5,5]]
dim        = 2

tb    = C.convertFile2PyTree('tbDouble.cgns')
tbox  = C.convertFile2PyTree('tboxDouble.cgns')
D_IBM._setDfar(tb,dFar)

t_AMR = G_AMR.generateAMRMesh(tb=tb, toffset=None, levelMax=levelMax, vmins=vmins, dim=dim,
                              check=False, opt=False, octreeMode=1, tbox=tbox, localDir=LOCAL)
test.testT(t_AMR,1)
#C.convertPyTree2File(t_AMR,'check_multi_t1_2D.cgns')
