# - integNorm (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

# Node mesh
m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
m = C.initVars(m, 'Density', 1.)
t = C.newPyTree(['Base',2]); t[2][1][2].append(m)
res = P.integNorm(t, 'Density'); print res
