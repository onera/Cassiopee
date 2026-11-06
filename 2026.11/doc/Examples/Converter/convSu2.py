# - convertArrays2File (SU2 format) -
import Converter as C

# Read mesh file
a = C.convertFile2Arrays("mesh_naca.su2", "fmt_su2")
C.convertArrays2File(a, "out.plt", "bin_tp")

# Rewrite mesh file
C.convertArrays2File(a, "out.su2", "fmt_su2")
