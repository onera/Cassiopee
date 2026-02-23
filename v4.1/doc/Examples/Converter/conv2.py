# - convertFile2Arrays (formatted tecplot) -
# - convertArrays2File (binary tecplot) -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("in360.tp", "fmt_tp")

# Write arrays to file
C.convertArrays2File(arrays, "out.plt", "bin_tp")
