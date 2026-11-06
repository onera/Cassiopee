# - transformTc2 (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.IBM as D_IBM
import Geom.PyTree as D
import numpy

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
a = C.node2Center(a)
for z in Internal.getZones(a):
    Internal._createChild(z, 'IBCD_2_'+z[0] , 'ZoneSubRegion_t', value=z[0])
    Internal._createChild(z, '2_IBCD_3_'+z[0] , 'ZoneSubRegion_t', value=z[0])

Nlength = numpy.zeros((10),numpy.float64)
for z in Internal.getZones(a):
    subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
    for zsr in subRegions:
        Internal._createChild(zsr, 'ZoneRole', 'DataArray_t', value='Donor')
        Internal._createChild(zsr, 'GridLocation', 'GridLocation_t', value='CellCenter')
        zsr[2].append(['Pressure', Nlength, [], 'DataArray_t'])
        zsr[2].append(['Density', Nlength, [], 'DataArray_t'])
        zsr[2].append(['VelocityX', Nlength, [], 'DataArray_t'])
        zsr[2].append(['VelocityY', Nlength, [], 'DataArray_t'])
        zsr[2].append(['VelocityZ', Nlength, [], 'DataArray_t'])

a=D_IBM.transformTc2(a)

C.convertPyTree2File(a, 'out.cgns')
