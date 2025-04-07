# - convertFile2Arrays (binary png) -
import Converter as C

a = C.convertFile2Arrays('Data/test.png', 'bin_png')
C.convertArrays2File(a, 'out.plt', 'bin_tp')
