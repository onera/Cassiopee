# - initBCDataSet (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
N = 10; NFACES = (N-1)*(N-1)
a = G.cart((0,0,0), (1,1,1), (N,N,N))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=NFACES*[1.], parent=d)

# Init all BCs
C._initBCDataSet(a, '{MomentumX}=2.*minimum({Density},0)')
# Init only BC node
bc = Internal.getNodeFromName(a, 'wall')
C._initBCDataSet(bc, 'MomentumY', 3.)
test.testT(a)
