# - convertArrays2File (GMSH file format) -
import Converter as C
import Generator as G

a = G.cartHexa( (0,0,0), (1,1,1), (3,3,3) )
b = G.cartTetra( (5,0,0), (1,1,1), (3,3,3) )

# write
C.convertArrays2File([a,b], 'out.msh', 'fmt_gmsh')

# Reread
a = C.convertFile2Arrays("out.msh", "fmt_gmsh")
