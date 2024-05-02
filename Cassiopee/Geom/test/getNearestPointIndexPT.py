# - getNearestPointIndex (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D

a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
inds = D.getNearestPointIndex(a, (0.55,0.34,0)); print(inds)
