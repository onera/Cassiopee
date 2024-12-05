# - convertFile2Arrays (fmt STL) -
import Converter as C
import Generator as G

a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = C.convertArrays2File([a], 'out.stl', 'fmt_stl')
a = C.convertFile2Arrays('out.stl', 'fmt_stl')
