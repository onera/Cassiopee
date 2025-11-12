# - addGordonSurface (array) -
import OCC
import numpy

N1 = 10
ucurves = []
for i in range(N1):
    hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
    points = numpy.array([[0,0,i],[1,1,i],[2,0,i]], dtype=numpy.float64, order='F')
    OCC.occ.addSpline(hook, points, 1)
    ucurves.append(hook)

N2 = 3
vcurves = []
hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
points = numpy.array([[0,0,0],[0,0,9]], dtype=numpy.float64, order='F')
OCC.occ.addSpline(hook, points, 1)
vcurves.append(hook)
hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
points = numpy.array([[1,1,0],[1,1,9]], dtype=numpy.float64, order='F')
OCC.occ.addSpline(hook, points, 1)
vcurves.append(hook)
hook = OCC.occ.createEmptyCAD("empty.stp", "fmt_step")
points = numpy.array([[2,0,0],[2,0,9]], dtype=numpy.float64, order='F')
OCC.occ.addSpline(hook, points, 1)
vcurves.append(hook)

hook = OCC.occ.mergeCAD(ucurves+vcurves)
OCC.occ.addGordonSurface(hook, [i for i in range(1,N1+1)], [i for i in range(N1+1,N1+N2+1)])

OCC.occ.writeCAD(hook, "out.step", "fmt_step")
