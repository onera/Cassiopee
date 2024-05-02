# - convertFile2Arrays (fmt v3d) -
# - convertArrays2File (fmt v3d) -
import Converter as C

# Lit le fichier in.plt dans arrays
arrays = C.convertFile2Arrays("infmt.v3d", "fmt_v3d")

# Sauvegarde la liste au format fmt_v3d
C.convertArrays2File(arrays, "outfmt.v3d", "fmt_v3d")
