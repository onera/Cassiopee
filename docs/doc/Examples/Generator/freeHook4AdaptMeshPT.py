# - freeHook4AdaptMesh (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

# Returns the hook on the tree structure
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,11))
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=3)
hooka = G._createHook4AdaptMesh(a, dim=3, splitInfos=None)
G.freeHook4AdaptMesh(hooka)
C.convertPyTree2File(a, "out.cgns")