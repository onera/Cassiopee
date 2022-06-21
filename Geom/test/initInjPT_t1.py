# - initInj (pyTree) -
import Geom.IBM as IBM
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy 
import KCore.test as test

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
a = C.node2Center(a)
for z in Internal.getZones(a):
    Internal._createChild(z, 'IBCD_5_'+z[0] , 'ZoneSubRegion_t', value=z[0])


Nlength = numpy.zeros((10),numpy.float64)
for z in Internal.getZones(a):
    subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
    for zsr in subRegions:
        Internal._createChild(zsr, 'ZoneRole', 'DataArray_t', value='Donor')
        Internal._createChild(zsr, 'GridLocation', 'GridLocation_t', value='CellCenter')


        zsr[2].append(['StagnationEnthalpy', Nlength, [], 'DataArray_t'])
        zsr[2].append(['StagnationPressure', Nlength, [], 'DataArray_t'])
        zsr[2].append(['dirx', Nlength, [], 'DataArray_t'])
        zsr[2].append(['diry', Nlength, [], 'DataArray_t'])
        zsr[2].append(['dirz', Nlength, [], 'DataArray_t'])

        
        Internal._createChild(zsr, 'FamilyName', 'FamilyName_t', value='CART_LOCAL')

a=IBM.initInj(a,'CART_LOCAL',10,20,injDir=[0.5,0.5,0.])
test.testT(a, 1)
