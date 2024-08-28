# - convertArrays2File (INRIA mesh format) -
import Converter as C

# Read mesh file
a = C.convertFile2Arrays("falcon.mesh", "fmt_mesh")
C.convertArrays2File(a, "out.plt", "bin_tp")

# Rewrite mesh file
C.convertArrays2File(a, "out.mesh", "fmt_mesh")
