# - convertFile2Arrays (binary 3DS) -
import Converter as C

a = C.convertFile2Arrays('Data/box.3ds', 'bin_3ds')
C.convertArrays2File(a, 'out.plt', 'bin_tp')
