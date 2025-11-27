# - getVolumeMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as T

# --- NGon grids ---
# method=1
a = G.cartNGon((0,0,0), (1,1,1), (3,5,7))
a = G.getVolumeMap(a, method=1)
C.convertPyTree2File(a, 'out.cgns')
T.testT(a, 1)

# method=0, api 1
a = G.cartNGon((0,0,0), (1,1,1), (3,5,7), api=1)
a = G.getVolumeMap(a, method=0)
vol = Internal.getNodeFromName(a, 'vol')[1]
T.testO(vol, 2)

# method=0, api 3
a = G.cartNGon((0,0,0), (1,1,1), (3,5,7), api=3)
a = G.getVolumeMap(a, method=0)
vol = Internal.getNodeFromName(a, 'vol')[1]
T.testO(vol, 2)
