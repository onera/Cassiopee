# - convertFile2Arrays (array) -
# Conversion d'un fichier tecplot 108 non-structure
# en fichier tecplot 75 non-structure
import Converter as C

a = C.convertFile2Arrays("Data/unstr3.plt", "bin_tp"); print(a)
C.convertArrays2File(a, "out.plt", "bin_tp")
