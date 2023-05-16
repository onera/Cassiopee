# - convertArrays2File (Foam format) -
import Converter.PyTree as C
import Generator.PyTree as G

t = G.cartNGon((0,0,0),(1,1,1),(21,21,2))
t = G.close(t)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall')
t = C.newPyTree(['Base', t])
#t = C.convertFile2PyTree("t.cgns")

# Write mesh file
C.convertPyTree2File(t, "out_foam", "fmt_foam")

# Reread
b = C.convertFile2PyTree("out_foam", "fmt_foam")
C.convertPyTree2File(b, 'foam_ngon.cgns')
