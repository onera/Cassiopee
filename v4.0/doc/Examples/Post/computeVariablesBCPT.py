# - computeVariablesBC (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
Internal.newDataArray('Density', value=10*[1.], parent=d)
Internal.newDataArray('MomentumX', value=10*[1.], parent=d)
Internal.newDataArray('MomentumY', value=10*[0.], parent=d)
Internal.newDataArray('MomentumZ', value=10*[0.], parent=d)
Internal.newDataArray('EnergyStagnationDensity', value=10*[1.], parent=d)

# Init all BCs
P._computeVariablesBC(a, ['Pressure'])

C.convertPyTree2File(a, 'out.cgns')
