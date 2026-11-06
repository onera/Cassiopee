# - getMassCenter (array) -
import Geom as D

a = D.sphere((1,0,0), 1., N=30)
Xc = D.getMassCenter(a)
print(Xc)