# - integ -
import Generator as G
import Converter as C
import Post as P
from numpy import *

nodes = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 11, 1))
C.convertArrays2File([nodes],'out.plt','bin_tp')

ni = nodes[2]
nj = nodes[3]

dens = ones( (1, ni * nj), float64 )
densa = ['t', dens, ni, nj, 1]

data = [ [nodes, densa ] ]

res = P.integNorm([nodes],[densa],[])
print(res)
