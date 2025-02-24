# - setInterpDataIBM 2D (pyTree) -
import Generator.PyTree as G
import Geom.IBM as DIBM
import Connector.IBM as XIBM
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
import Transform.PyTree as T
import Converter.Mpi as Cmpi

# --- Front Type=1
N = 50; h = 1./(N-1)
c = G.cylinder((0.,0.,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); c[0] = 'cyl'
DIBM._setIBCType(c, 'Musker')
tb = C.newPyTree(['CYL', Internal.getZones(c)])

a = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=2)
XIBM._blankingIBM(t, tb, dimPb=2)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = XIBM.buildFrontIBM(t, tc, dimPb=2, cartesian=False)
XIBM._setInterpDataIBM(t, tc, tb, front, dimPb=2, frontType=1, cartesian=False)
####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
Internal._rmNodesByName(tc, '.Solver#Param')
####
test.testT(t, 1)
test.testT(tc, 11)

# --- Front Type=42
N = 50; h = 1./(N-1)
c = G.cylinder((0.,0.,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); c[0] = 'cyl'
DIBM._setIBCType(c, 'Musker')
tb = C.newPyTree(['CYL', Internal.getZones(c)])

a = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=2, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=2, frontType=42)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = XIBM.buildFrontIBM(t, tc, dimPb=2, frontType=42, cartesian=False)
XIBM._setInterpDataIBM(t, tc, tb, front, dimPb=2, frontType=42, cartesian=False)
####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
Internal._rmNodesByName(tc, '.Solver#Param')
####
test.testT(t, 2)
test.testT(tc, 21)

# --- Front Type=42, correction multicorps=True
N = 50; h = 1./(N-1)
c1 = G.cylinder((-0.5,-0.5,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); c1[0] = 'cyl1'
c2 = G.cylinder((0.5,0.5,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); c2[0] = 'cyl2'
tb = C.newPyTree(['CYL1', Internal.getZones(c1),
                  'CYL2', Internal.getZones(c2)])
DIBM._setIBCType(tb, 'Musker')

a = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=2, correctionMultiCorpsF42=True, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=2, correctionMultiCorpsF42=True, frontType=42)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = XIBM.buildFrontIBM(t, tc, dimPb=2, frontType=42, cartesian=False)
XIBM._setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=2, frontType=42, cartesian=False)
####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
Internal._rmNodesByName(tc, '.Solver#Param')
####
test.testT(t, 3)
test.testT(tc, 31)

# --- Front Type=42 , blankingF42=True
N = 50; h = 1./(N-1)
c = G.cylinder((0.,0.,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); c[0] = 'cyl'
DIBM._setIBCType(c, 'Musker')
tb = C.newPyTree(['CYL', Internal.getZones(c)])

a = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=2, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=2, frontType=42, blankingF42=True)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = XIBM.buildFrontIBM(t, tc, dimPb=2, cartesian=False)
XIBM._setInterpDataIBM(t, tc, tb, front, dimPb=2, frontType=42, cartesian=False)
####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
Internal._rmNodesByName(tc, '.Solver#Param')
####
test.testT(t, 4)
test.testT(tc, 41)

# --- Front Type=42, twoFronts=True
N = 50; h = 1./(N-1)
c = G.cylinder((0.,0.,0.), 0.2, 1., 360., 0., 10., (2*N,1,1)); c[0] = 'cyl'
DIBM._setIBCType(c, 'Musker')
tb = C.newPyTree(['CYL', Internal.getZones(c)])

a = G.cart((-1,-1,0), (h,h,1), (2*N,2*N,1)); a[0] = 'cart'
t = C.newPyTree(['CART', Internal.getZones(a)])
T._addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
C._fillEmptyBCWith(t, 'farfield', 'BCFarfield')

XIBM._dist2wallIBM(t, tb, dimPb=2, frontType=42)
XIBM._blankingIBM(t, tb, dimPb=2, frontType=42, twoFronts=True)
tc = C.node2Center(t)
t, tc, front, front2, frontWMM = XIBM.buildFrontIBM(t, tc, dimPb=2, frontType=42, twoFronts=True, cartesian=False)
XIBM._setInterpDataIBM(t, tc, tb, front, front2=front2, dimPb=2, frontType=42, twoFronts=True, cartesian=False)
####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
Internal._rmNodesByName(tc, '.Solver#Param')
####
test.testT(t, 5)
test.testT(tc, 51)
