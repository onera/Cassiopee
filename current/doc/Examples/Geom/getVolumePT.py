# - getVolume (pyTree) -
import Geom.PyTree as D
import math

a = D.sphere((0,0,0), 1., N=30)
vol = D.getVolume(a)
print(vol, 4./3.*math.pi*1*1*1)
