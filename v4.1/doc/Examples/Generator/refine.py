# - refine (array) -
import Generator as G
import Converter as C

a = G.cart( (0,0,0), (0.1,0.1,0.1), (20,20,1) )

a = G.refine(a, 1.5, 1)
C.convertArrays2File([a], 'out.plt')
