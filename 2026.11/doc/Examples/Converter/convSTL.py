# - convertFile2Arrays (binary STL) -
import Converter as C

a = C.convertFile2Arrays('Data/piece-a-percer.stl', 'bin_stl')
C.convertArrays2File(a, 'out.plt', 'bin_tp')
