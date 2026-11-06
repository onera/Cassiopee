# - createHook4AdaptMesh (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Returns the hook on the tree structure
a = G.cartHexa((0,0,0),(0.1,0.1,0.1),(11,11,11))
C._fillEmptyBCWith(a, 'nref', 'BCFarfield', dim=3)
C._initVars(a, "F", 1.)
C._initVars(a, '{centers:indicator}=({centers:CoordinateX})>0.5')
hooka = G._createHook4AdaptMesh(a, dim=3, splitInfos=None)
G.freeHook4AdaptMesh(hooka)
test.testT(a, 1)