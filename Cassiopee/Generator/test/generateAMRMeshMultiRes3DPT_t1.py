import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.IBM as D_IBM
import Geom.PyTree as D
import Generator.AMR as G_AMR
import KCore.test as test

LOCAL      = test.getLocal()
snear      = 0.2
dFar       = 5
levelMax   = 2
vmins      = [[6,6]]
vminsTbox  = [[5,5]]
dim        = 3

a = D.sphere((0,0,0), 0.5, 20)
b = D.cylinder((1.5,0,-0.5), 0.5, 1., N=10, ntype='QUAD')
tb = C.newPyTree(["BODYSPHERE",a, 'BODYCYL',b])
D_IBM._setSnear(tb,snear)
D_IBM._setDfar(tb,dFar)
D_IBM._setIBCType(tb,"Musker")
#C.convertPyTree2File(tb, 'check_tb.cgns')

tbox= None

t_AMR = G_AMR.generateAMRMesh(tb=tb, toffset=None, levelMax=levelMax, vmins=vmins, dim=dim,
                              check=True, opt=False, octreeMode=1, tbox=tbox, localDir=LOCAL)
test.testT(t_AMR,1)
#C.convertPyTree2File(t_AMR,'check_multi_t1_3D.cgns')
