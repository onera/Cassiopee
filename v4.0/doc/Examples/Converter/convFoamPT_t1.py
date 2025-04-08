# - convertPyTree2File (Foam format) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

t = G.cartNGon((0,0,0),(1,1,1),(21,21,2))
t = G.close(t)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall')
t = C.newPyTree(['Base', t])

# Write mesh file
C.convertPyTree2File(t, LOCAL+"/out_foam", "fmt_foam")

# Reread
b = C.convertFile2PyTree(LOCAL+"/out_foam", "fmt_foam")

#C.convertPyTree2File(b, 'out.cgns')
test.testT(b, 1)
