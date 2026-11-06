# - extractIBMWallSolution (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test
import Post.IBM as P_IBM
import Geom.IBM as D_IBM
import Connector.IBM as X_IBM
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Mpi as Cmpi
import numpy

def _updateFlowSolutionIBM(tc):
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:4] == "IBCD":
                XPW = Internal.getNodeFromName1(zsr, 'CoordinateX_PW')[1]
                YPW = Internal.getNodeFromName1(zsr, 'CoordinateY_PW')[1]
                ZPW = Internal.getNodeFromName1(zsr, 'CoordinateZ_PW')[1]

                Internal.getNodeFromName(zsr, 'Pressure')[1] = 2*XPW + 3*YPW + 4*ZPW
                Internal.getNodeFromName(zsr, 'Density')[1] = numpy.ones_like(XPW)
                Internal.getNodeFromName(zsr, 'VelocityX')[1] = 2*numpy.ones_like(XPW)
                Internal.getNodeFromName(zsr, 'VelocityY')[1] = 3*numpy.ones_like(XPW)
                Internal.getNodeFromName(zsr, 'VelocityZ')[1] = 4*numpy.ones_like(XPW)
    return None

def splitAndDispatchTree(t):
    t = T.splitNParts(t, Cmpi.size)
    for i, z in enumerate(Internal.getZones(t)):
        Cmpi._setProc(z, i)
    Cmpi._convert2PartialTree(t, Cmpi.rank)

    return t

#### 2D (Circle)

# Create t_in
N = 51
z = G.cart((-2,-2, 0),(4./(N-1),4./(N-1),1),(N,N,1))
t = C.newPyTree(['Base', z])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
t = splitAndDispatchTree(t)

# Create case tree
z = D.circle((0.,0.,0.), 1., N=101)
tb = C.newPyTree(['Base', z])
D_IBM._setSnear(tb, 0.05)
D_IBM._setDfar(tb, 10.)
D_IBM._setIBCType(tb, 'Musker')
C._addState(tb, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=2, GoverningEquations='NSTurbulent')
C._addState(t, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=2, GoverningEquations='NSTurbulent')

# prep IBM + init IBM info in tc
X_IBM._dist2wallIBM(t, tb, dimPb=2)
X_IBM._blankingIBM(t, tb, dimPb=2, frontType=1, IBCType=1, depth=2)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = X_IBM.buildFrontIBM(t, tc, tb=tb, dimPb=2, frontType=1)
X_IBM._setInterpDataIBM(t, tc, tb, front=front, dimPb=2, frontType=1, IBCType=1, depth=2)

_updateFlowSolutionIBM(tc)

# test 1: all
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=2, ibctypes=[], famZones=[])
P_IBM._computeSkinVariables(ts, tc, graphIBCDPost, dimPb=2, ibctypes=[], famZones=[])
if Cmpi.master: test.testT(ts, 1)

#### 3D (Sphere)

# Create t_in
N = 51
z = G.cart((-10,-10,-10),(20./(N-1),20./(N-1),20./(N-1)),(N,N,N))
t = C.newPyTree(['Base', z])
t = splitAndDispatchTree(t)

# Create case tree
z = D.sphere6((0.,0.,0.), 1., N=25, ntype='TRI')
tb = C.newPyTree(['Base', z])
D_IBM._setSnear(tb, 0.1)
D_IBM._setDfar(tb, 10.)
D_IBM._setIBCType(tb, 'Musker')
C._addState(tb, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=3, GoverningEquations='NSTurbulent')
C._addState(t, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=3, GoverningEquations='NSTurbulent')

# prep IBM + init IBM info in tc
X_IBM._dist2wallIBM(t, tb, dimPb=3)
X_IBM._blankingIBM(t, tb, dimPb=3, frontType=1, IBCType=1, depth=2)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = X_IBM.buildFrontIBM(t, tc, tb=tb, dimPb=3, frontType=1)
X_IBM._setInterpDataIBM(t, tc, tb, front=front, dimPb=3, frontType=1, IBCType=1, depth=2)
_updateFlowSolutionIBM(tc)

# test 2: all
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=3, ibctypes=[], famZones=[])
P_IBM._computeSkinVariables(ts, tc, graphIBCDPost, dimPb=3, ibctypes=[], famZones=[])
if Cmpi.master: test.testT(ts, 2)