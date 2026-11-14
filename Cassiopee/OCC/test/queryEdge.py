# - queryEdge (array) -
import OCC

# line
hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (1,0,0))
d = OCC.queryEdge(hook, 1)
print(d)

# circle
hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 2.)
d = OCC.queryEdge(hook, 1)
print(d)

# spline1
hook = OCC.createEmptyCAD()
import Geom
a = Geom.polyline([(0,0,0), (1,1,0), (2,0,0)])
OCC._addSpline(hook, a[1], method=0, degree=3)
d = OCC.queryEdge(hook, 1)
print(d)

# ellipse
hook = OCC.createEmptyCAD()
OCC._addEllipse(hook, (0,0,0), (0,0,1), 2., 1.)
d = OCC.queryEdge(hook, 1)
print(d)
