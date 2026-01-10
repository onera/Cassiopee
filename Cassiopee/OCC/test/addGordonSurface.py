# - addGordonSurface (array) -
import OCC
import numpy

N1 = 10
ucurves = []
for i in range(N1):
    hook = OCC.createEmptyCAD()
    points = numpy.array([[0,0,i],[1,1,i],[2,0,i]], dtype=numpy.float64, order='F')
    OCC._addSpline(hook, points, 1, 3)
    ucurves.append(hook)

N2 = 3
vcurves = []
hook = OCC.createEmptyCAD()
points = numpy.array([[0,0,0],[0,0,9]], dtype=numpy.float64, order='F')
OCC._addSpline(hook, points, 1, 3)
vcurves.append(hook)
hook = OCC.createEmptyCAD()
points = numpy.array([[1,1,0],[1,1,9]], dtype=numpy.float64, order='F')
OCC._addSpline(hook, points, 1, 3)
vcurves.append(hook)
hook = OCC.createEmptyCAD()
points = numpy.array([[2,0,0],[2,0,9]], dtype=numpy.float64, order='F')
OCC._addSpline(hook, points, 1, 3)
vcurves.append(hook)

hook = OCC.mergeCAD(ucurves+vcurves)
OCC._addGordonSurface(hook, [i for i in range(1,N1+1)], [i for i in range(N1+1,N1+N2+1)])

OCC.writeCAD(hook, "out.step")
