# - join (array) -
import Transform as T
import Converter as C
import Generator as G

a1 = G.cartTetra((0.,0.,0.), (1.,1.,1), (11,11,1))
a2 = G.cartTetra((10.,0.,0.), (1.,1.,1), (10,10,1))
a = T.join(a1, a2)
C.convertArrays2File([a], 'out.plt')
