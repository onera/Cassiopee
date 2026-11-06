# - attractEdges (array) -
import Converter as C
import Transform as T
import Generator as G

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1) )
b = G.cartTetra( (-1,0,0), (0.5, 1, 1), (30,1,1) )
b = T.rotate(b, (0,0,0), (0,0,1), 30.)

C.convertArrays2File([a,b], 'out.plt')
