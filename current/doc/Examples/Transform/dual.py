# - dual (arrays) -
import Converter as C
import Generator as G
import Transform as T

ni = 5; nj = 5; nk = 1
a = G.cart((0,0,0),(1,1,1),(ni,nj,nk))
a = C.convertArray2NGon(a); a = G.close(a)
res = T.dual(a)
C.convertArrays2File([res],'out.tp')
