# - rotate (array) -
import Generator as G
import Transform as T
import Converter as C
import math
import KCore.test as test
ni = 91; nj = 21
a = G.cylinder((0,0,0),0.5,1.,0.,90.,1.,(ni,nj,1))
a = C.initVars(a,'vx', 0.)
a = C.initVars(a,'vy', 0.)
a = C.initVars(a,'vz', 1.)
for j in range(nj):
    for i in range(ni):
        i0 = i * math.pi/2./90.
        a[1][3,i+j*ni] = math.cos(i0)
        a[1][4,i+j*ni] = math.sin(i0)

# Rotate with an axis and an angle
b = T.rotate(a, (0.,0.,0.), (0.,0.,1.), 90., vectors=[['vx','vy','vz']])
test.testA([b],1)
# Rotate with axis transformations
c = T.rotate(a, (0.,0.,0.), ((1.,0.,0.),(0,1,0),(0,0,1)),
             ((1,1,0), (1,-1,0), (0,0,1)), vectors=[['vx','vy','vz']] )
test.testA([c],2)

# Rotate with three angles
d = T.rotate(a, (0.,0.,0.), (0.,0.,90.), vectors=[['vx','vy','vz']])
test.testA([d],3)
