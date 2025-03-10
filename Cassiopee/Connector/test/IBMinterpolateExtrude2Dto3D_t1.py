# extrude 2D mesh to 3D with Cartesian approach & prepareIBMDataExtrude
import Apps.Fast.IBM as App
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

#C.convertPyTree2File(t3D ,LOCAL+'/t3D_checking.cgns')
#C.convertPyTree2File(tc3D,LOCAL+'/tc3D_checking.cgns')

####Keep below  - however, certain bug fixes in App.prepare1 are required
#### 1) recompute sign distance needs to be corrected (Nk[z[0]] += 2*ific-1 +c  # -1)
#### 2) tb extrude needs to be corrected (if not self.isFilamentOnly: X_IBM._signDistance(t))
##
#### ============================================
#### Previous Extrude & Interpolation 3D Approach
#### ============================================
##
##X._applyBCOverlaps(t2, depth=2, loc='centers', val=1, cellNName='cellN')
##C._initVars(t2, '{centers:TurbulentDistanceAllBC}={centers:TurbulentDistance}')
##C._initVars(t2, '{centers:cellNIBC_blank}=0*({centers:cellN}<1)+1*({centers:cellN}>0)')
##C._initVars(t2, '{centers:cellNIBC_hole}={centers:cellN}')
##X._applyBCOverlaps(t2, depth=2, loc='centers', val=2, cellNName='cellN')
##
#### Extrusion for 3D Mesh
##t3Dorig, tb3Dorig = App.extrudeCartesian(t2, bodySurface2, extrusion=extrusion, NPas=NPas, span=span, dz=span/(NPas-1), isAutoPeriodic=True)
##for t in [t3Dorig,tb3Dorig]:
##    zmax   = C.getMaxValue(t, 'CoordinateZ');
##    zmin   = C.getMinValue(t, 'CoordinateZ');
##    zavg   = (zmax+zmin)/2
##    T._translate(t, (0,0,0-zavg))
##
##BCs = Internal.getNodesFromType(t3Dorig, "BC_t")
##for bc in BCs:
##    if Internal.getValue(bc)=='BCautoperiod':
##        nameSplit=bc[0].split(".")
##        lenLocal=len(nameSplit)
##        bc[0]='.'.join(nameSplit[0:lenLocal-1])
##
### Interpolation 3D
##interpDataType    = 0
##order             = 2
##t3Dorig, tc3Dorig = App.prepare1(tb3Dorig, None, None, t_in=t3Dorig, extrusion=extrusion, interpDataType=interpDataType, order=order)
##
##C._rmVars(t3Dorig, ['centers:TurbulentDistanceAllBC','centers:cellNIBC_blank','centers:cellNIBC_hole'])
##Internal._rmNodesByName(t3Dorig, 'ReferenceState')
##Internal._rmNodesByName(t3Dorig, 'FlowEquationSet')
##
##Internal._rmNodesByName(tc3Dorig, 'ReferenceState')
##Internal._rmNodesByName(tc3Dorig, 'FlowEquationSet')
##
##C._rmVars(tc3Dorig, ['TurbulentDistanceAllBC','cellNIBC_blank','cellNIBC_hole'])
##
##test.testT(t3Dorig  ,1)
##test.testT(tc3Dorig ,2)
##
###C.convertPyTree2File(t3Dorig , LOCAL+'/t3Dorig_checking.cgns')
###C.convertPyTree2File(tc3Dorig, LOCAL+'/tc3Dorig_checking.cgns')
