# - rotate (array) -
import Generator as G
import Transform as T
import Converter as C
import numpy
import KCore.test as test

ni = 91; nj = 21
a = G.cylinder((0,0,0),0.5,1.,0.,90.,1.,(ni,nj,1))
a = C.initVars(a,'vx', 0.)
a = C.initVars(a,'vy', 0.)
a = C.initVars(a,'vz', 1.)

i0 = numpy.arange(ni) * numpy.pi / 180.
a[1][3,:] = numpy.tile(numpy.cos(i0), nj)
a[1][4,:] = numpy.tile(numpy.sin(i0), nj)

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
