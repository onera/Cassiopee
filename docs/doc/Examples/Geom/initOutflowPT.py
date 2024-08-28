# - initOutflow (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.IBM as D_IBM
import Geom.PyTree as D
import numpy

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
a = C.node2Center(a)
for z in Internal.getZones(a):
    Internal._createChild(z, 'IBCD_4_'+z[0] , 'ZoneSubRegion_t', value=z[0])


Nlength = numpy.zeros((10),numpy.float64)
for z in Internal.getZones(a):
    subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
    for zsr in subRegions:
        Internal._createChild(zsr, 'ZoneRole', 'DataArray_t', value='Donor')
        Internal._createChild(zsr, 'GridLocation', 'GridLocation_t', value='CellCenter')
        zsr[2].append(['Pressure', Nlength, [], 'DataArray_t'])
        Internal._createChild(zsr, 'FamilyName', 'FamilyName_t', value='CART_LOCAL')

a=D_IBM.initOutflow(a,'CART_LOCAL',101325)

C.convertPyTree2File(a, 'out.cgns')
