import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.IBM as D_IBM
import Generator.AMR as G_AMR
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

LOCAL      = test.getLocal()
snear      = 0.001*5
dFar       = 10
levelMax   = 20
vmins      = [[15,14,14,10,11,5,4]]
vminsTbox  = [[5,5]]
dim        = 2

# NACAs
a = D.naca(12., N=1001)
dimPb = 2
D_IBM._setSnear(a, snear)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 10.)
b = T.translate(a, (0.,-0.70,0.))
D_IBM._setSnear(b, snear*4)
tb = C.newPyTree(["BODYNACA",a, 'COPYNACA', b])
D_IBM._setDfar(tb,dFar)

# Tboxs
l1 = D.line((1.55,-0.20, 0.0), (2.30,-0.20, 0.0), 200)
l2 = D.line((2.30,-0.20, 0.0), (2.30, 0.85, 0.0), 200)
l3 = D.line((2.30, 0.85, 0.0), (1.55, 0.85, 0.0), 200)
l4 = D.line((1.55, 0.85, 0.0), (1.55,-0.20, 0.0), 200)
a  = T.join([l1,l2,l3,l4])
D_IBM._setSnear(a, snear*4)

b  = T.translate(a, (0.,-1.30,0.))
D_IBM._setSnear(b, snear*8)
tbox = C.newPyTree(["BODY1",a, 'COPY1', b])
D_IBM._setDfar(tbox,20)
D_IBM._setIBCType(tbox,"None")

t_AMR = G_AMR.generateAMRMesh(tb=tb, toffset=None, levelMax=levelMax, vmins=vmins, dim=dim,
                              check=False, opt=False, octreeMode=1, tbox=tbox, localDir=LOCAL)
test.testT(t_AMR,1)
#C.convertPyTree2File(t_AMR,'check_multi_t1_2D.cgns')
