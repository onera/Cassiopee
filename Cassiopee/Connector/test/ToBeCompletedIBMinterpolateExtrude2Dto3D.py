
## DO NOT RUN THIS TEST CASE - IT IS IN PROGRESS AND WILL BE COMPLETED SHORTLY

# extrude 2D mesh to 3D with Cartesian approach
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Converter.Internal as Internal
import Generator.IBM as G_IBM
import Geom.IBM as Geom
import Connector.IBM as X_IBM
import KCore.test as test
import sys

LOCAL           = test.getLocal()
bodySurfaceFile = LOCAL  + '/naca0012.cgns'

Lcharac     = 0.03362355
Lz          = 0.2*Lcharac
dfar        = 10.*Lcharac
vmin        = 16;
snears      = 1;

###Generate 2D Mesh
t2D,tc2D=X_IBM.prepareIBMDataPara(bodySurfaceFile     , None                  , None           ,
                                  snears=snears       , dfar=dfar             , vmin=vmin      )
test.testT(t2D ,1)
test.testT(tc2D,2)
#C.convertPyTree2File(t2D,LOCAL+'/t2D_checking.cgns')

####Extrusion for 3D Mesh
bodySurface = C.convertFile2PyTree(bodySurfaceFile)

extrusion   = 'cart'
span        = Lz
NPas        = 10+1 #number of nodes
t3D, tb3D   = G_IBM.extrudeCartesianZDir(t2D, bodySurface, extrusion=extrusion, NPas=NPas, span=span, dz=span/(NPas-1), isAutoPeriodic=True)

for t in [t3D,tb3D]:
    zmax   = C.getMaxValue(t, 'CoordinateZ');
    zmin   = C.getMinValue(t, 'CoordinateZ');
    zavg   = (zmax+zmin)/2
    T._translate(t, (0,0,0-zavg))
test.testT(t3D  ,3)
test.testT(tb3D ,4)

#C.convertPyTree2File(t3D,LOCAL+'/t3D_checking.cgns')
#C.convertPyTree2File(tb3D,LOCAL+'/tb3D_checking.cgns')

#####Interpolation 3D
t3D, tc3D      = X_IBM.prepareIBMDataExtrude(tb3D, None, None, t3D, extrusion=extrusion)
test.testT(t3D  ,5)
test.testT(tc3D ,6)

