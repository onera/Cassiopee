# extrude 2D mesh to 3D with Cartesian approach & prepareIBMDataExtrude
#import Apps.Fast.IBM as App
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Generator.IBM as G_IBM
import Connector.IBM as X_IBM
import Connector.PyTree as X
import KCore.test as test
import sys, os
import numpy


LOCAL           = test.getLocal()
bodySurfaceFile = 'naca1DNS.cgns'

# Prepare
vmin      = 42
dfars     = 5
snears    = 1
t, tc = X_IBM.prepareIBMData(bodySurfaceFile, None         , None     ,
                             snears=snears  , dfars=dfars  , vmin=vmin,
                             check=False    , frontType=1  , cartesian=False)

#C.convertPyTree2File(t,LOCAL+'/t2D_checking.cgns')
#C.convertPyTree2File(tc,LOCAL+'/tc2D_checking.cgns')

bodySurface = C.convertFile2PyTree(bodySurfaceFile)

t2           = Internal.copyTree(t)
bodySurface2 = Internal.copyTree(bodySurface)

## ===========================================
## Current Extrude & Interpolation 3D Approach
## ===========================================

# Extrusion for 3D Mesh
extrusion   = 'cart'
span        = 1
NPas        = 4+1 #number of nodes
t3D, tb3D   = G_IBM.extrudeCartesianZDir(t, bodySurface, extrusion=extrusion, NPas=NPas, span=span, dz=span/(NPas-1), isAutoPeriodic=True)

for t in [t3D,tb3D]:
    zmax   = C.getMaxValue(t, 'CoordinateZ');
    zmin   = C.getMinValue(t, 'CoordinateZ');
    zavg   = (zmax+zmin)/2
    T._translate(t, (0,0,0-zavg))

#C.convertPyTree2File(t3D ,LOCAL+'/t3D_checking.cgns')
#C.convertPyTree2File(tb3D,LOCAL+'/tb3D_checking.cgns')

# Interpolation 3D
t3D, tc3D      = X_IBM.prepareIBMDataExtrude(tb3D, None, None, t3D, extrusion=extrusion)

####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
for b in Internal.getBases(tc3D):
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

Internal._rmNodesByName(t3D, 'FlowEquationSet')
Internal._rmNodesByName(t3D, 'ReferenceState')

Internal._rmNodesByName(tc3D, 'FlowEquationSet')
Internal._rmNodesByName(tc3D, 'ReferenceState')

test.testT(t3D  ,1)
test.testT(tc3D ,2)