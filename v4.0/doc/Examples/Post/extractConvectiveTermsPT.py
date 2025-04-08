# - extractConvectiveTerms (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Post.IBM as P_IBM
import copy
import numpy

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
a = C.node2Center(a)
for z in Internal.getZones(a):
    Internal._createChild(z, 'IBCD_2_'+z[0] , 'ZoneSubRegion_t', value=z[0])

Nlength  = numpy.zeros((10),numpy.float64)
for z in Internal.getZones(a):
    subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
    for zsr in subRegions:
        Internal._createChild(zsr, 'ZoneRole', 'DataArray_t', value='Donor')
        Internal._createChild(zsr, 'GridLocation', 'GridLocation_t', value='CellCenter')

        zsr[2].append(['CoordinateX_PW', copy.copy(Nlength), [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PW', copy.copy(Nlength), [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PW', copy.copy(Nlength), [], 'DataArray_t'])

        zsr[2].append(['CoordinateX_PC', copy.copy(Nlength)+14, [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PC', copy.copy(Nlength)+14, [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PC', copy.copy(Nlength)+14, [], 'DataArray_t'])

        zsr[2].append(['CoordinateX_PI', copy.copy(Nlength)+15, [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PI', copy.copy(Nlength)+15, [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PI', copy.copy(Nlength)+15, [], 'DataArray_t'])

        zsr[2].append(['Pressure', Nlength+3, [], 'DataArray_t'])
        zsr[2].append(['Density' , copy.copy(Nlength)+1, [], 'DataArray_t'])

        zsr[2].append(['gradxPressure', copy.copy(Nlength)+2, [], 'DataArray_t'])
        zsr[2].append(['gradyPressure', copy.copy(Nlength)+2, [], 'DataArray_t'])
        zsr[2].append(['gradzPressure', copy.copy(Nlength)+2, [], 'DataArray_t'])

        zsr[2].append(['gradxVelocityX', copy.copy(Nlength)+3, [], 'DataArray_t'])
        zsr[2].append(['gradyVelocityX', copy.copy(Nlength)+3, [], 'DataArray_t'])
        zsr[2].append(['gradzVelocityX', copy.copy(Nlength)+3, [], 'DataArray_t'])

        zsr[2].append(['gradxVelocityY', copy.copy(Nlength)+4, [], 'DataArray_t'])
        zsr[2].append(['gradyVelocityY', copy.copy(Nlength)+4, [], 'DataArray_t'])
        zsr[2].append(['gradzVelocityY', copy.copy(Nlength)+4, [], 'DataArray_t'])

        zsr[2].append(['gradxVelocityZ', copy.copy(Nlength)+5, [], 'DataArray_t'])
        zsr[2].append(['gradyVelocityZ', copy.copy(Nlength)+5, [], 'DataArray_t'])
        zsr[2].append(['gradzVelocityZ', copy.copy(Nlength)+5, [], 'DataArray_t'])

        zsr[2].append(['VelocityX', copy.copy(Nlength)+6, [], 'DataArray_t'])
        zsr[2].append(['VelocityY', copy.copy(Nlength)+7, [], 'DataArray_t'])
        zsr[2].append(['VelocityZ', copy.copy(Nlength)+8, [], 'DataArray_t'])

a=P_IBM.extractConvectiveTerms(a)
C.convertPyTree2File(a,'out.cgns')
