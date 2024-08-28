# - collapse (array) -
import Converter as C
import Generator as G
import Transform as T

a = G.cartTetra((0.,0.,0.),(0.1,0.01,1.),(20,2,1))
b = T.collapse(a)
C.convertArrays2File([a,b], "out.plt")
