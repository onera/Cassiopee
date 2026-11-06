# - convertArrays2File (Cedre fmt format) -
import Converter as C
import Generator as G

# Create grid
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
C.convertArrays2File([a], "out.d", "fmt_cedre")

# Reread file
a = C.convertFile2Arrays("out.d", "fmt_cedre")
