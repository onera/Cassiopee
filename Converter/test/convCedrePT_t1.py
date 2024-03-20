# - convertPyTree2File (Cedre fmt format) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Create grid
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
C._fillEmptyBCWith(a, 'wall', 'BCWall')

# write file
C.convertPyTree2File(a, LOCAL+"/out.d", "fmt_cedre")

# Reread file
a = C.convertFile2PyTree(LOCAL+"/out.d", "fmt_cedre")

#C.convertPyTree2File(a, 'out.cgns')
test.testT(a, 1)
