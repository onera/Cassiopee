# extrude 2D mesh to 3D with Cartesian approach
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Generator.IBM as G_IBM
import Connector.IBM as X_IBM
import KCore.test as test
import sys

LOCAL           = test.getLocal()
bodySurfaceFile = '../../Connector/test/naca1DNS.cgns'

# Prepare
vmin      = 42
dfars     = 5
snears    = 1
t, tc = X_IBM.prepareIBMData(bodySurfaceFile, None         , None     ,
                             snears=snears  , dfars=dfars  , vmin=vmin,
                             check=False    , frontType=1  , cartesian=False)

####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
for b in Internal.getBases(tc):
    for z in Internal.getZones(b):
        pos = 0
        z2 = Internal.copyRef(z)
        for zs in z2[2]:
            if 'ID' in zs[0] or 'IBCD' in zs[0]:
                Internal.addChild(z, zs, pos)
                pos +=2
            else:
                pos += 1
####

test.testT(t ,1)
test.testT(tc,2)
#C.convertPyTree2File(t,LOCAL+'/t2D_checking.cgns')


####Extrusion for 3D Mesh
bodySurface = C.convertFile2PyTree(bodySurfaceFile)

t2           = Internal.copyTree(t)
bodySurface2 = Internal.copyTree(bodySurface)

extrusion   = 'cart'
span        = 1
NPas        = 10+1 #number of nodes
t3D, tb3D   = G_IBM.extrudeCartesianZDir(t, bodySurface, extrusion=extrusion, NPas=NPas, span=span, dz=span/(NPas-1), isAutoPeriodic=True)

for t in [t3D,tb3D]:
    zmax   = C.getMaxValue(t, 'CoordinateZ');
    zmin   = C.getMinValue(t, 'CoordinateZ');
    zavg   = (zmax+zmin)/2
    T._translate(t, (0,0,0-zavg))
C._rmVars(t3D, ['centers:cellN'])
test.testT(t3D  ,3)
test.testT(tb3D ,4)

#C.convertPyTree2File(t3D ,LOCAL+'/t3D_checking.cgns')
#C.convertPyTree2File(tb3D,LOCAL+'/tb3D_checking.cgns')
#

t3Dorig, tb3Dorig = App.extrudeCartesian(t2, bodySurface2, extrusion=extrusion, NPas=NPas, span=span, dz=span/(NPas-1), isAutoPeriodic=True)
for t in [t3Dorig,tb3Dorig]:
    zmax   = C.getMaxValue(t, 'CoordinateZ');
    zmin   = C.getMinValue(t, 'CoordinateZ');
    zavg   = (zmax+zmin)/2
    T._translate(t, (0,0,0-zavg))

BCs = Internal.getNodesFromType(t3Dorig, "BC_t")
for bc in BCs:
    if Internal.getValue(bc)=='BCautoperiod':
        nameSplit=bc[0].split(".")
        lenLocal=len(nameSplit)
        bc[0]='.'.join(nameSplit[0:lenLocal-1])
C._rmVars(t3Dorig, ['centers:cellN'])
test.testT(t3Dorig,3)

#C.convertPyTree2File(t3Dorig ,LOCAL+'/t3Dorig_checking.cgns')
