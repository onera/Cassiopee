# - splitNParts (array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0,0,0), (1,1,1), (101,101,41))
b = G.cart((10,0,0), (1,1,1), (121,61,81))
c = G.cart((20,0,0), (1,1,1), (101,61,131))
res = T.splitNParts([a,b,c], 32, multigrid=0, dirs=[1,2,3])

C.convertArrays2File(res, 'out.plt')
