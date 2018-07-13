# - convertFile2Arrays (binaire v3d) -
# - convertArrays2File (formatted and binary v3d) -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("inbin.v3d", "bin_v3d")

# Write a binary file with endian conversion
C.convertArrays2File(arrays, "out.v3d", "fmt_v3d")
