# - fillHole (array) -
import OCC

# method 0 (planar sketch)
hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)
OCC._fillHole(hook, [1], [], 0)
OCC.writeCAD(hook, "out.step")

# method 1 (general sketch - circle)
hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)
OCC._fillHole(hook, [1], [], 1)
OCC.writeCAD(hook, "out2.step")

# method 1 (general sketch - square)
hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (10,0,0))
OCC._addLine(hook, (10,0,0), (10,10,3))
OCC._addLine(hook, (10,10,3), (0,10,0))
OCC._addLine(hook, (0,10,0), (0,0,0))
OCC._fillHole(hook, [1,2,3,4], [], 1)
OCC.writeCAD(hook, "out3.step")

# method 1 (general sketch - with spline)
hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (10,0,0))
import Geom
points = Geom.polyline([(10,0,0), (5,5,4), (10,10,3)])
OCC._addSpline(hook, points[1], method=0, degree=3)
OCC._addLine(hook, (10,10,3), (0,10,0))
OCC._addLine(hook, (0,10,0), (0,0,0))
OCC._fillHole(hook, [1,2,3,4], [], 1)
OCC.writeCAD(hook, "out4.step")
