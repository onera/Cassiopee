# - convertFile2Arrays (tecplot binaire 108) -
# - convertArrays2File (plot3d) -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("in360.plt", "bin_tp")

# Write a binary file
C.convertArrays2File(arrays, "out.dat.gbin,out.dat.qbin",
                     "bin_plot3d",['int',4])
