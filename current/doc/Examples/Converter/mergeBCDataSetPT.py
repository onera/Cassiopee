# - mergeBCDataSet (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G

N = 10; NFACES = (N-1)*(N-1)
a = G.cart((0,0,0), (1,1,1), (N,N,N))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=NFACES*[1.], parent=d)
d2 = Internal.newBCDataSet(name='BCDataSet2', value='UserDefined',
                           gridLocation='FaceCenter', parent=b)
d2 = Internal.newBCData('BCNeumann', parent=d2)
d2 = Internal.newDataArray('VelocityX', value=NFACES*[0.3], parent=d2)
C._mergeBCDataSets(a)
C.convertPyTree2File(a, 'out.cgns')
