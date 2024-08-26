import Converter.PyTree as C
import Post.IBM as P_IBM
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.Internal as Internal
import numpy
import KCore.test as test

##### 2D TESTS ############################################################
a = D.circle((0,0,0),1,N=40); a=C.convertArray2Tetra(a); a = G.close(a)
tb = C.newPyTree(["Base",a])

#C.convertPyTree2File(tb, 'tb.cgns')

x0 = -2; N = 41; h = 2*abs(x0)/(N-1)
z = G.cart((x0,x0,x0),(h,h,h),(N,N,1))
zname = Internal.getName(z)
zsr = Internal.createNode('IBCD'+'_3_'+zname, 'ZoneSubRegion_t', value=zname)
Internal._createChild(zsr, 'GridLocation', 'GridLocation_t', value='CellCenter')
# mimic the IBM wall pt info
a2 = D.circle((0,0,0),1, N=40); a2 = C.convertArray2Tetra(a2); a2 = G.close(a2)
GC = Internal.getNodeFromType(a2,"GridCoordinates_t")
FSN = Internal.getNodeFromType(a2,'FlowSolution_t')
nIBC = Internal.getZoneDim(a2)[1]
XP = numpy.zeros((nIBC),numpy.float64)
XN = Internal.getNodeFromName(GC,'CoordinateX')[1]; XP[:]=XN[:]
YP = numpy.zeros((nIBC),numpy.float64)
YN = Internal.getNodeFromName(GC,'CoordinateY')[1]; YP[:]=YN[:]
ZP = numpy.ones((nIBC),numpy.float64)*0.005
DENS = numpy.ones((nIBC),numpy.float64)
PRESS = 101325*numpy.ones((nIBC),numpy.float64)
zsr[2].append(['CoordinateX_PW', XP, [], 'DataArray_t'])
zsr[2].append(['CoordinateY_PW', YP, [], 'DataArray_t'])
zsr[2].append(['CoordinateZ_PW', ZP, [], 'DataArray_t'])
zsr[2].append(['Pressure', PRESS, [], 'DataArray_t'])
zsr[2].append(['Density', DENS, [], 'DataArray_t'])
z[2].append(zsr)
tc = C.newPyTree(['CART']); tc[2][1][2].append(z)

#C.convertPyTree2File(tc, 'tc.cgns')

graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=2, prepareMLS=False)
test.testT(ts, 1)

graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=2, prepareMLS=True)
test.testT(ts, 2)

##### 3D TESTS ############################################################

a = D.sphere((0,0,0),1,N=30); a=C.convertArray2Tetra(a); a = G.close(a)
tb = C.newPyTree(["Base",a])

x0 = -2; N = 41; h = 2*abs(x0)/(N-1)
z = G.cart((x0,x0,x0),(h,h,h),(N,N,N))
zname = Internal.getName(z)
zsr = Internal.createNode('IBCD'+'_3_'+zname, 'ZoneSubRegion_t', value=zname)
Internal._createChild(zsr, 'GridLocation', 'GridLocation_t', value='CellCenter')

# mimic the IBM wall pt info
a2 = D.sphere((0,0,0),1, N=30); a2 = C.convertArray2Tetra(a2); a2 = G.close(a2)
GC = Internal.getNodeFromType(a2,"GridCoordinates_t")
FSN = Internal.getNodeFromType(a2,'FlowSolution_t')
nIBC = Internal.getZoneDim(a2)[1]
XP = numpy.zeros((nIBC),numpy.float64)
XN = Internal.getNodeFromName(GC,'CoordinateX')[1]; XP[:]=XN[:]
YP = numpy.zeros((nIBC),numpy.float64)
YN = Internal.getNodeFromName(GC,'CoordinateY')[1]; YP[:]=YN[:]
ZP = numpy.zeros((nIBC),numpy.float64)
ZN = Internal.getNodeFromName(GC,'CoordinateZ')[1]; ZP[:]=ZN[:]
DENS = numpy.ones((nIBC),numpy.float64)
PRESS = 101325*numpy.ones((nIBC),numpy.float64)
zsr[2].append(['CoordinateX_PW', XP, [], 'DataArray_t'])
zsr[2].append(['CoordinateY_PW', YP, [], 'DataArray_t'])
zsr[2].append(['CoordinateZ_PW', ZP, [], 'DataArray_t'])
zsr[2].append(['Pressure', PRESS, [], 'DataArray_t'])
zsr[2].append(['Density', DENS, [], 'DataArray_t'])
z[2].append(zsr)
tc = C.newPyTree(['CART']); tc[2][1][2].append(z)

graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=3, prepareMLS=False)
test.testT(ts, 3)

graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=3, prepareMLS=True)
test.testT(ts, 4)
