# - convertFile2Arrays (binary v3d) -
# - convertArrays2File (formatted and binary v3d) -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("infmt.v3d", "fmt_v3d")

# Write a binary file with endian conversion
C.convertArrays2File(arrays, "out.v3d", "bin_v3d", endian='big')

# Write a binary file with "elsA v3d format"
C.convertArrays2File(arrays, "out_elsa.v3d", "fmt_v3d", dataFormat='%14.7e')
