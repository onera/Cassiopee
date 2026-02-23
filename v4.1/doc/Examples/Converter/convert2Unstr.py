# - convertFile2Arrays (bin_tp unstructured) -
import Converter as C

# lit un fichier non structure
# le met dans un array
# et le remet dans un fichier non structure
uarrays = C.convertFile2Arrays("Data/unstr.plt", "bin_tp")

arrays = C.convertFile2Arrays("dauphin.plt", "bin_tp")
arrays.append(uarrays[0])
C.convertArrays2File(arrays, "out.plt", "bin_tp")
