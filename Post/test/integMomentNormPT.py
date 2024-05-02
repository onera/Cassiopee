# - integMomentNorm (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
m = C.initVars(m, 'Density', 1.)
res = P.integMomentNorm(m, var='Density', center=(5.,5.,0.)); print(res)
#>> [[-3.685500000000002, 3.6855, 0.0]]
