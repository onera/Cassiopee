# - convertFile2Arrays (formatted Obj) -
import Converter as C

a = C.convertFile2Arrays('Data/cube.obj', 'fmt_obj')
C.convertArrays2File(a, 'out.plt', 'bin_tp')
