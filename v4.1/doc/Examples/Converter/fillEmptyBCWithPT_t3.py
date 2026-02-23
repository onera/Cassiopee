# - fillEmptyBC (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Transform.PyTree as T
import KCore.test as test

# BE - ajoute en tant que Element_t
a = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
C._fillEmptyBCWith(a, 'wallv', 'BCWallViscous')
test.testT(a, 1)

# BE + initial BC
#a = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))
#p = P.exteriorFaces(a)
#p = T.splitSharpEdges(p)
#C._addBC2Zone(a, 'wall', 'BCWall', subzone=p[0])
#C._fillEmptyBCWith(a, 'wallv', 'BCWallViscous')
#test.testT(a, 2)

# ME to complete
