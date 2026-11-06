# - convertFile2Arrays (binary plot3d) -
# - convertArrays2File (binary tecplot) -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("out.dat", "bin_plot3d")

# Write arrays to file
C.convertArrays2File(arrays, "out.plt", "bin_tp")
