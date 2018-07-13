# - collapse (array) -
import Converter as C
import Generator as G
import Transform as T

hi = 0.1; hj = 0.01; hk = 1.
ni = 20; nj = 2; nk = 1
a = G.cartTetra((0.,0.,0.),(hi,hj,hk),(ni,nj,nk))
b = T.collapse(a)
C.convertArrays2File([a,b], "out.plt")
