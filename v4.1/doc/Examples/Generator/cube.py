# - cube -
import Converter as C
import Generator.Shapes as S

m = S.square((0,0), (1,1), 0.5, (10,10,5))

# Ecriture fichier
C.convertArrays2File(m, "out.plt", "bin_tp")
