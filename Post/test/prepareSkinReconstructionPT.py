# - prepareSkinReconstruction (pyTree) -
import Converter.PyTree as C
import Post.IBM as P_IBM
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.Internal as Internal
import numpy

a = D.sphere((0,0,0),1,N=20); a=C.convertArray2Tetra(a); a = G.close(a)
ts = C.newPyTree(["Base",a])

x0 = -2; N = 41; h = 2*abs(x0)/(N-1)
z = G.cart((x0,x0,x0),(h,h,h),(N,N,N))
zname = Internal.getName(z)
zsr = Internal.createNode('IBCD_'+zname, 'ZoneSubRegion_t', value=zname)
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

tl, graphs = P_IBM._prepareSkinReconstruction(ts,tc)
C.convertPyTree2File(tl,"out.cgns")
