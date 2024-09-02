# - perturbate (array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = T.perturbate(a, 0.1)
C.convertArrays2File(a, "out.plt")
