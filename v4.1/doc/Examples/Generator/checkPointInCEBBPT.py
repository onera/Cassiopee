# - checkPointInCEBB (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T

Ni = 20; Nj = 20
a2 = G.cart((-0.1,0,0),(0.5/Ni, 0.5/Nj, 1), (Ni,Nj,2))
a2 = T.rotate(a2, (-0.1,0,0), (0,0,1), 0.22)

# Check if point is in CEBB of a2
val = G.checkPointInCEBB(a2, (0.04839, 0.03873, 0.5)); print(val)
