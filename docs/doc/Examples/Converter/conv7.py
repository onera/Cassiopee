# - convertFile2Arrays : bin_tp non structure -
# - convertArrays2File : bin_tp non structure -
import Converter as C

unstArrays = C.convertFile2Arrays("Data/unstr.plt", "bin_tp")
C.convertArrays2File(unstArrays, "out.plt", "bin_tp")
