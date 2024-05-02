# - convertFile2Arrays (binary tecplot) -
# - convertArrays2File (formatted tecplot) -
import Converter as C

# Lit le fichier in.plt dans arrays
arrays = C.convertFile2Arrays("in.plt", "bin_tp")

# Sauvegarde la liste au format fmt_tp
C.convertArrays2File(arrays, "out.tp", "fmt_tp")
