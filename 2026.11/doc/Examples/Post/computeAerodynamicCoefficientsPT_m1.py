# - computeAerodynamicCoefficients (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Geom.PyTree as D
import Post.IBM as P_IBM
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import KCore.test as test
import numpy

def _initSurfaceFlowSolution(tb, case=0):
    if case == 0: # no friction
        for z in Internal.getZones(tb):
            C._initVars(z, '{Pressure} = 3*{CoordinateX} + 6*{CoordinateY} + 9*{CoordinateZ} + 20')
            C._initVars(z, '{Density} = 1.')
            C._initVars(z, '{VelocityX} = 0.')
            C._initVars(z, '{VelocityY} = 0.')
            C._initVars(z, '{VelocityZ} = 0.')
            C._initVars(z, '{utau} = 0.')
            C._initVars(z, '{gradxPressure} = 0.')
            C._initVars(z, '{gradyPressure} = 0.')
            C._initVars(z, '{gradzPressure} = 0.')
    else:
        for z in Internal.getZones(tb):
            C._initVars(z, '{Pressure} = 3*{CoordinateX} + 6*{CoordinateY} + 9*{CoordinateZ} + 20')
            C._initVars(z, '{Density} = 3*{CoordinateX} + 6*{CoordinateY} + 9*{CoordinateZ} + 20')
            C._initVars(z, '{VelocityX} = -{CoordinateY}')
            C._initVars(z, '{VelocityY} = {CoordinateX}')
            C._initVars(z, '{VelocityZ} = 0.')
            C._initVars(z, '{utau} = 1.')
            C._initVars(z, '{gradxPressure} = 0.')
            C._initVars(z, '{gradyPressure} = 0.')
            C._initVars(z, '{gradzPressure} = 0.')
    return None

def _initReferenceState(tb):
    C._addState(tb, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=2, GoverningEquations='NSTurbulent')

    refstate = Internal.getNodeByType2(tb, 'ReferenceState_t')

    Internal.getNodeFromName1(refstate, 'Pressure')[1] = 1.
    Internal.getNodeFromName1(refstate, 'Density')[1] = 2.

    Vx = Internal.getNodeFromName1(refstate, 'VelocityX')[1]
    Vy = Internal.getNodeFromName1(refstate, 'VelocityY')[1]
    Vz = Internal.getNodeFromName1(refstate, 'VelocityZ')[1]

    Vnorm = numpy.sqrt(Vx**2 + Vy**2 + Vy**2)

    Internal.getNodeFromName1(refstate, 'VelocityX')[1] = Vx/Vnorm
    Internal.getNodeFromName1(refstate, 'VelocityY')[1] = Vy/Vnorm
    Internal.getNodeFromName1(refstate, 'VelocityZ')[1] = Vz/Vnorm

    return None

def splitAndDispatchTree(t):
    t = T.splitNParts(t, Cmpi.size)
    for i, z in enumerate(Internal.getZones(t)):
        Cmpi._setProc(z, i)
    Cmpi._convert2PartialTree(t, Cmpi.rank)

    return t

#### 2D (Circle)

z = D.circle((0.,0.,0.), 1., N=500)
tb = C.newPyTree(['Base', z])

tb = splitAndDispatchTree(tb)
_initSurfaceFlowSolution(tb, case=1)
_initReferenceState(tb)

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=None, dimPb=2, famZones=[], Pref=None, center=(0.,0.,0.), verbose=1)
if Cmpi.master:
    test.testT(wall, 11)
    test.testO(aeroLoads, 12)
# -> Fx_P, Fy_P, Fz_P = 3*pi, 6*pi, 0
# -> Fx_V, Fy_V, Fz_V = -6*pi, 3*pi, 0
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 40*pi
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=2, Sref=None, Lref=None, Qref=None, alpha=0., beta=0., verbose=1)
if Cmpi.master:
    test.testT(wall, 21)
    test.testO(aeroLoads, 22)
# -> Fx_P, Fy_P, Fz_P = 1.5, 0, 3
# -> Fx_V, Fy_V, Fz_V = -3, 0, 1.5
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 20, 0

#### 3D (Sphere)

z = D.sphere6((0.,0.,0.), 1., N=25, ntype='TRI')
tb = C.newPyTree(['Base', z])

tb = splitAndDispatchTree(tb)
_initSurfaceFlowSolution(tb, case=1)
_initReferenceState(tb)

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=None, dimPb=3, famZones=[], Pref=None, center=(0.,0.,0.), verbose=1)
if Cmpi.master:
    test.testT(wall, 31)
    test.testO(aeroLoads, 32)
# -> Fx_P, Fy_P, Fz_P = -4*pi, -8*pi, -12*pi
# -> Fx_V, Fy_V, Fz_V = -3*pi**2, (3*pi**2)/2, 0
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 20*pi**2
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=3, Sref=None, Lref=None, Qref=None, alpha=0., beta=0., verbose=1)
if Cmpi.master:
    test.testT(wall, 41)
    test.testO(aeroLoads, 42)
# -> Fx_P, Fy_P, Fz_P = -1, -2, -3
# -> Fx_V, Fy_V, Fz_V = -3*pi/4, 3*pi/8, 0
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 5*pi