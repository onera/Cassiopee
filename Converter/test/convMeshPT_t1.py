# - convertPyTree2File (mesh fmt format) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Create grid
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
#C._fillEmptyBCWith(a, 'wall', 'BCWall')

# write file
C.convertPyTree2File(a, "out.mesh", "fmt_mesh")

# Reread file
a = C.convertFile2PyTree("out.mesh", "fmt_mesh")

#C.convertPyTree2File(a, 'out.cgns')
test.testT(a, 1)
