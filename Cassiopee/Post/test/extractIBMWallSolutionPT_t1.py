# - extractIBMWallSolution (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test
import Post.IBM as P_IBM
import Geom.IBM as D_IBM
import Connector.IBM as X_IBM
import Post.PyTree as P
import Generator.PyTree as G
import Transform.PyTree as T
import numpy

def _updateFlowSolutionIBM(tc):
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:4] == "IBCD":
                Internal._rmNodesByName(zsr, 'utau')
                Internal._rmNodesByName(zsr, 'yplus')

                XPW = Internal.getNodeFromName1(zsr, 'CoordinateX_PW')[1]
                YPW = Internal.getNodeFromName1(zsr, 'CoordinateY_PW')[1]
                ZPW = Internal.getNodeFromName1(zsr, 'CoordinateZ_PW')[1]

                Internal.getNodeFromName(zsr, 'Pressure')[1] = 2*XPW + 3*YPW + 4*ZPW
                Internal.getNodeFromName(zsr, 'Density')[1] = numpy.ones_like(XPW)
                Internal.getNodeFromName(zsr, 'VelocityX')[1] = 2*numpy.ones_like(XPW)
                Internal.getNodeFromName(zsr, 'VelocityY')[1] = 3*numpy.ones_like(XPW)
                Internal.getNodeFromName(zsr, 'VelocityZ')[1] = 4*numpy.ones_like(XPW)


    return None

#### 2D (NACA0012)

# Create t_in
Nx = 151
Ny = 51
z = G.cart((-1,-0.5, 0),(3./(Nx-1),1./(Ny-1),1),(Nx,Ny,1))
t = C.newPyTree(['Base', z])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)

# Create naca zone
z = D.naca('0012', N=101)

z_le = P.selectCells(z, '{CoordinateX} < 0.1', strict=1)
z_le[0] = 'naca_le' # leading edge
C._tagWithFamily(z_le, 'wing')
D_IBM._setIBCType(z_le, 'Musker')

z_main = P.selectCells(z, '{CoordinateX} >= 0.1')
z_main[0] = 'naca_main' # main wing
C._tagWithFamily(z_main, 'wing')
D_IBM._setIBCType(z_main, 'slip')

# Create circle zone
z2 = D.circle((-0.15,0.,0.), 0.05, N=101)
z2[0] = 'circle' # circle obstacle
C._tagWithFamily(z2, 'obstacle')
D_IBM._setIBCType(z2, 'Musker')

# Create case tree
tb = C.newPyTree(['NACA', [z_le, z_main], 'CIRCLE', [z2]])
D_IBM._setSnear(tb, 0.05)
D_IBM._setDfar(tb, 10.)
C._addState(tb, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=2, GoverningEquations='NSTurbulent')
C._addState(t, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=2, GoverningEquations='NSTurbulent')

# prep IBM + init IBM info in tc
a, ac = X_IBM.prepareIBMData(tb, t_in=t, t_out=None, tc_out=None, frontType=1, check=False)
_updateFlowSolutionIBM(ac)

# test 1: All + prepareMLS = False
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=2, ibctypes=[], famZones=[], prepareMLS=False)
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=2, ibctypes=[], famZones=[])
test.testT(ts, 1)

# test 2: All + prepareMLS = True
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=2, ibctypes=[], famZones=[])
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=2, ibctypes=[], famZones=[])
test.testT(ts, 2)

# test 3: Musker IBC only
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=2, ibctypes=[3], famZones=[])
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=2, ibctypes=[3], famZones=[])
test.testT(ts, 3)

# test 4: naca family zones only
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=2, ibctypes=[], famZones=['wing'])
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=2, ibctypes=[], famZones=['wing'])
test.testT(ts, 4)

# test 5: circle family zones only
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=2, ibctypes=[], famZones=['obstacle'])
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=2, ibctypes=[], famZones=['obstacle'])
test.testT(ts, 5)

# test 6: Musker IBC + naca family zones only
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=2, ibctypes=[3], famZones=['wing'])
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=2, ibctypes=[3], famZones=['wing'])
test.testT(ts, 6)

#### 3D (SPHERE)

# Create t_in
N = 51
z = G.cart((-10,-10,-10),(20./(N-1),20./(N-1),20./(N-1)),(N,N,N))
t = C.newPyTree(['Base', z])

# Create case tree
z = D.sphere6((0.,0.,0.), 1., N=25, ntype='TRI')
tb = C.newPyTree(['Base', z])
D_IBM._setSnear(tb, 0.1)
D_IBM._setDfar(tb, 10.)
D_IBM._setIBCType(tb, 'Musker')
C._addState(tb, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=3, GoverningEquations='NSTurbulent')
C._addState(t, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6, EquationDimension=3, GoverningEquations='NSTurbulent')

# prep IBM + init IBM info in tc
a, ac = X_IBM.prepareIBMData(tb, t_in=t, t_out=None, tc_out=None, frontType=1, check=False)
_updateFlowSolutionIBM(ac)

# test 7: All + prepare MLS = False
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=3, ibctypes=[], famZones=[], prepareMLS=False)
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=3, ibctypes=[], famZones=[])
test.testT(ts, 7)

# test 8: All + prepare MLS = True
graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, ac, dimPb=3, ibctypes=[], famZones=[], prepareMLS=True)
P_IBM._computeSkinVariables(ts, ac, graphIBCDPost, dimPb=3, ibctypes=[], famZones=[])
test.testT(ts, 8)