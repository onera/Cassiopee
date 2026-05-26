# - distance (array) -
import Geom as D

a = D.circle((0,0,0), 1., N=30)
b = D.circle((0,0,0), 1.1, N=30)
d = D.distance(a, b)
print(d)

a = D.sphere((0,0,0), 1., N=30)
b = D.sphere((0,0,0), 1.1, N=30)
d = D.distance(a, b)
print(d)
