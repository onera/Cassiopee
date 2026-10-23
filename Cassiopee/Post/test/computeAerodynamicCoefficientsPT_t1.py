# - computeAerodynamicCoefficients (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Geom.PyTree as D
import Post.IBM as P_IBM
import Post.PyTree as P
import Generator.PyTree as G
import Transform.PyTree as T
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

#### 2D (Rectangle)
print("----- CASE 1 = RECTANGLE -----")

z = G.cartHexa((0,0,0), (0.1,0.1,0.), (41, 11, 1))
z = P.exteriorFaces(z)
z = D.uniformize(z, 500)
tb = C.newPyTree(['Base', z])

_initSurfaceFlowSolution(tb, case=0)
_initReferenceState(tb)

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=None, dimPb=2, famZones=[], Pref=1, center=(1.,0.5,0.), verbose=1)
# -> Fx_P, Fy_P, Fz_P = 12, 24, 0
# -> Fx_V, Fy_V, Fz_V = 0, 0, 0 
# -> Mx_P, My_P, Mz_P = 0, 0, 24
# -> Mx_V, My_V, Mz_V = 0, 0, 0
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=2, Sref=4, Lref=2, Qref=1, alpha=0., beta=0., verbose=1)
# -> Fx_P, Fy_P, Fz_P = 3, 0, 6
# -> Fx_V, Fy_V, Fz_V = 0, 0, 0 
# -> Mx_P, My_P, Mz_P = 0, 3, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 0

#### 2D (Rotated rectangle - angle = 15deg)
print("\n----- CASE 2 = ROTATED RECTANGLE -----")

z = G.cartHexa((0,0,0), (0.1,0.1,0.), (41, 11, 1))
z = P.exteriorFaces(z)
z = D.uniformize(z, 500)
tb = C.newPyTree(['Base', z])

_initSurfaceFlowSolution(tb, case=0)
_initReferenceState(tb)

T._rotate(tb, (1,0.5,0), (0,0,1), 15.)

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=None, dimPb=2, famZones=[], Pref=1, center=(1.,0.5,0.), verbose=1)
# -> Fx_P, Fy_P, Fz_P = 12*cos(15)-24*sin(15), 0, 12*sin(15)+24*cos(15)
# -> Fx_V, Fy_V, Fz_V = 0, 0, 0 
# -> Mx_P, My_P, Mz_P = 0, 0, 24
# -> Mx_V, My_V, Mz_V = 0, 0, 0
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=2, Sref=1, Lref=1, Qref=1, alpha=15, beta=0., verbose=1)
# -> Fx_P, Fy_P, Fz_P = 12, 0, 24
# -> Fx_V, Fy_V, Fz_V = 0, 0, 0 
# -> Mx_P, My_P, Mz_P = 0, 24, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 0

#### 2D (Circle)
print("\n----- CASE 3 = CIRCLE -----")

z = D.circle((0.,0.,0.), 1., N=500)
tb = C.newPyTree(['Base', z])

_initSurfaceFlowSolution(tb, case=1)
_initReferenceState(tb)

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=None, dimPb=2, famZones=[], Pref=None, center=(0.,0.,0.), verbose=1)
# -> Fx_P, Fy_P, Fz_P = 3*pi, 6*pi, 0
# -> Fx_V, Fy_V, Fz_V = -6*pi, 3*pi, 0 
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 40*pi
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=2, Sref=None, Lref=None, Qref=None, alpha=0., beta=0., verbose=1)
# -> Fx_P, Fy_P, Fz_P = 1.5, 0, 3
# -> Fx_V, Fy_V, Fz_V = -3, 0, 1.5
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 20, 0

#### 2D (Circle + 2nd tree to extract the pressure solution)
print("\n----- CASE 4 = CIRCLE WITH 2ND TREE FOR PRESSURE -----")

z = D.circle((0.,0.,0.), 1., N=500)
tb = C.newPyTree(['Base', z])

_initSurfaceFlowSolution(tb, case=1)
_initReferenceState(tb)

tb2 = Internal.copyTree(tb) # 2nd circle - used to extract the pressure solution
for z2 in Internal.getZones(tb2):
    C._initVars(z2, '{Pressure} = 6*{CoordinateX} + 12*{CoordinateY} + 18*{CoordinateZ} + 40')
    C._initVars(z2, '{gradxPressure} = 1.')
    C._initVars(z2, '{gradyPressure} = 2.')
    C._initVars(z2, '{gradzPressure} = 3.')

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=tb2, dimPb=2, famZones=[], Pref=None, center=(0.,0.,0.), verbose=1)
# -> Fx_P, Fy_P, Fz_P = 6*pi, 12*pi, 0
# -> Fx_V, Fy_V, Fz_V = -6*pi, 3*pi, 0 
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 40*pi
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=2, Sref=None, Lref=None, Qref=None, alpha=0., beta=0., verbose=1)
# -> Fx_P, Fy_P, Fz_P = 3, 0, 6
# -> Fx_V, Fy_V, Fz_V = -3, 0, 1.5
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 20, 0

#### 3D (Sphere)
print("\n----- CASE 5 = SPHERE -----")

z = D.sphere6((0.,0.,0.), 1., N=25, ntype='TRI')
tb = C.newPyTree(['Base', z])

_initSurfaceFlowSolution(tb, case=1)
_initReferenceState(tb)

wall, aeroLoads = P_IBM.computeAerodynamicLoads(tb, ts2=None, dimPb=3, famZones=[], Pref=None, center=(0.,0.,0.), verbose=1)
# -> Fx_P, Fy_P, Fz_P = -4*pi, -8*pi, -12*pi
# -> Fx_V, Fy_V, Fz_V = -3*pi**2, (3*pi**2)/2, 0
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 20*pi**2
wall, aeroLoads = P_IBM.computeAerodynamicCoefficients(wall, aeroLoads, dimPb=3, Sref=None, Lref=None, Qref=None, alpha=0., beta=0., verbose=1)
# -> Fx_P, Fy_P, Fz_P = -1, -2, -3
# -> Fx_V, Fy_V, Fz_V = -3*pi/4, (3*pi**2)/8, 0
# -> Mx_P, My_P, Mz_P = 0, 0, 0
# -> Mx_V, My_V, Mz_V = 0, 0, 5*pi