# - loft (array) -
import OCC
import numpy
from math import sin, cos

N = 5
ucurves = []
for i in range(N):
    hook = OCC.createEmptyCAD()
    OCC._addCircle(hook, (20,0,0), (0,0,1), 1.)
    OCC._rotate(hook, (0,0,0), (0,1,0), i/(N-1)*90.)
    ucurves.append(hook)
hook = OCC.mergeCAD(ucurves)

# loft without guides
OCC._loft(hook, [i for i in range(1,N+1)], [])
OCC.writeCAD(hook, "out.step")

# loft with guides
N2 = 2
vcurves = []
hook = OCC.createEmptyCAD()
R = 19.; a = numpy.pi/2.*1./4.
x1 = R*sin(0.); z1 = -R*cos(0.)
x2 = R*sin(a); z2 = -R*cos(a)
x3 = R*sin(2*a); z3 = -R*cos(2*a)
x4 = R*sin(3*a); z4 = -R*cos(3*a)
x5 = R*sin(4*a); z5 = -R*cos(4*a)
points = numpy.array([[x1,0,z1],[x2,0,z2],[x3,0,z3],[x4,0,z4],[x5,0,z5]], dtype=numpy.float64, order='F')
OCC._addSpline(hook, points, 1, 3)
vcurves.append(hook)
hook = OCC.createEmptyCAD()
R = 21.; a = numpy.pi/2.*1./4.
x1 = R*sin(0.); z1 = -R*cos(0.)
x2 = R*sin(a); z2 = -R*cos(a)
x3 = R*sin(2*a); z3 = -R*cos(2*a)
x4 = R*sin(3*a); z4 = -R*cos(3*a)
x5 = R*sin(4*a); z5 = -R*cos(4*a)
points = numpy.array([[x1,0,z1],[x2,0,z2],[x3,0,z3],[x4,0,z4],[x5,0,z5]], dtype=numpy.float64, order='F')
OCC._addSpline(hook, points, 1, 3)
vcurves.append(hook)

hook = OCC.mergeCAD(ucurves+vcurves)
OCC._loft(hook, [i for i in range(1,N+1)], [i for i in range(N+1, N+N2+1)])
OCC.writeCAD(hook, "out.step")
