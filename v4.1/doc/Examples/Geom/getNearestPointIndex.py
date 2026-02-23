# - getNearestPointIndex (array) -
import Generator as G
import Geom as D

a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
inds = D.getNearestPointIndex(a, (0.55,0.34,0)); print(inds)
inds = D.getNearestPointIndex(a, [(0.55,0.34,0), (0.56,0.32,0)]); print(inds)
