# - convertPyTree2File (SU2 fmt format) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

# Create grid
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
C._fillEmptyBCWith(a, 'wall', 'BCWall')

# write file
C.convertPyTree2File(a, LOCAL+"/out.su2", "fmt_su2")

# Reread file - pour l'instant les BCs ne sont pas relues
a = C.convertFile2PyTree(LOCAL+"/out.su2", "fmt_su2")

#C.convertPyTree2File(a, 'out.cgns')
test.testT(a, 1)
