# - exteriorFacesStructured (array) -
import Converter as C
import Post as P
import Generator as G

a = G.cart((0,0,0), (1,1,1), (4,4,6))
A = P.exteriorFacesStructured(a)
C.convertArrays2File(A, 'out.plt')
