# - getVolume (array) -
import Geom as D
import numpy

a = D.sphere((0,0,0), 1., N=30)
vol = D.getVolume(a)
print(vol, 4./3.*numpy.pi*1*1*1)
