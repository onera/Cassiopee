# - sweep (array) -
import OCC
import numpy
from math import sin, cos

# profile
hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (20,0,0), (0,0,1), 1.)

# path
R = 19.; a = numpy.pi/2.*1./4.
x1 = R*sin(0.); z1 = -R*cos(0.)
x2 = R*sin(a); z2 = -R*cos(a)
x3 = R*sin(2*a); z3 = -R*cos(2*a)
x4 = R*sin(3*a); z4 = -R*cos(3*a)
x5 = R*sin(4*a); z5 = -R*cos(4*a)
points = numpy.array([[x1,0,z1],[x2,0,z2],[x3,0,z3],[x4,0,z4],[x5,0,z5]], dtype=numpy.float64, order='F')
OCC._addSpline(hook, points, 1, 3)

OCC._sweep(hook, [1], [2])
OCC.writeCAD(hook, "out.step")
