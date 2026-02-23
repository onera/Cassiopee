# - convertArrays2File (fmt_pov) -
import Converter as C

a = C.convertFile2Arrays("dauphin_skin.plt", "bin_tp")

# fmt_pov only supports triangle meshes
b = C.convertArray2Tetra(a)

# This converts coordinates and Density field to fmt_pov
# using colormap 1 (c1). c0 means no color pigment.
C.convertArrays2File(b, "out.pov", "fmt_pov", ['colormap',1])

# Reread this file
a = C.convertFile2Arrays("out.pov", "fmt_pov")
C.convertArrays2File(a, "out.plt", "bin_tp")

# Execute povray
import os; os.system("povray -W800 -H600 +a0.3 +SP16 render.pov +P")
