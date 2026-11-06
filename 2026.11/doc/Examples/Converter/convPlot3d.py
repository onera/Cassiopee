# - convertFile2Arrays (binary plot3d) -
# - convertArrays2File (binary plot3d) -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("in.dat", "bin_plot3d")
arrays = C.addVars(arrays, ["ro", "rou", "rov", "row", "roE"])

# Write arrays to file
C.convertArrays2File(arrays, "out.dat", "bin_plot3d", ['int',8])
C.convertArrays2File(arrays, "out.plt", "bin_tp")
