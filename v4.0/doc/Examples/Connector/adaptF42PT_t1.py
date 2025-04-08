# - adaptF42 (prepareIBMData) (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.PyTree as D
import Geom.IBM as D_IBM
import Connector.IBM as X_IBM
import Generator.IBMmodelHeight as G_IBM_Height
import KCore.test as test
import numpy

def _initYplusTargetPoints(tc):
    hmod = G_IBM_Height.computeModelisationHeight(6.e6, Cf_law='ANSYS', yplus=2000., L=1.)

    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:4] == 'IBCD':
                yplus = Internal.getNodeFromName(zsr, 'yplus')
                if yplus is not None:
                    yplus = yplus[1]

                    XPW = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                    YPW = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                    ZPW = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                    XPC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                    YPC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                    ZPC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                    distCW = numpy.sqrt( (XPW-XPC)*(XPW-XPC) + (YPW-YPC)*(YPW-YPC) + (ZPW-ZPC)*(ZPW-ZPC))

                    Internal.getNodeFromName(zsr, 'yplus')[1] = distCW*(4000/hmod)
    return None

body = D.naca(12.,N=101)
tb  = C.newPyTree(['Base', body])
C._addState(tb, 'EquationDimension', 2)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')
C._addState(tb, adim='adim1', MInf=0.15, alphaZ=0., alphaY=0., ReInf=6.e6)
D_IBM._setIBCType(tb, 'Musker')
D_IBM._setSnear(tb, 0.001)
D_IBM._setDfar(tb, 25.)
t,tc = X_IBM.prepareIBMData(tb, None, None, check=False, frontType=42, yplus=2000, vmin=21, cartesian=False, expand=4)

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

test.testT(tc,1)
C.__ZoneNameServer__ = {}; C.__BCNameServer__ = {}; C.__BaseNameServer__ = {}
_initYplusTargetPoints(tc)
wallAdapt = X_IBM.createWallAdapt(tc)
t,tc = X_IBM.prepareIBMData(tb, None, None, check=False, frontType=42, yplus=2000, vmin=21, wallAdaptF42=wallAdapt, cartesian=False, expand=4)

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

test.testT(tc,2)
