# - integMoment (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

m = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
C._initVars(m,'vx',1.); C._initVars(m,'vy',0.); C._initVars(m,'vz',0.)
res = P.integMoment(m, center=(5.,5., 0.), vector=['vx','vy','vz']); print(res)
#>> [0.0, 0.0, 3.6854999999999976]
