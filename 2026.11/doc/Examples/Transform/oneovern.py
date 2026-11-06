# - oneovern (array) -
import Transform as T
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,1))
a2 = T.oneovern(a, (2,2,2))
C.convertArrays2File(a2, "out.plt")
