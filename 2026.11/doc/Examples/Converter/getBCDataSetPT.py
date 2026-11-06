# - getBCDataSet (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=10*[1.], parent=d)

# Get data array node list
datas = Internal.getBCDataSet(a, b); print(datas)
