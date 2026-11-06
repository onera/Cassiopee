# - convertFile2Arrays (binaire v3d) -
# - convertArrays2UsurpGeneric -
import Converter as C

# Read a file into arrays
arrays = C.convertFile2Arrays("dauphin_skin.plt", "bin_tp")

C.convertArrays2UsurpGeneric(arrays, arrays)
