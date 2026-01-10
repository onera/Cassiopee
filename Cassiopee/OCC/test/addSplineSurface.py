# - addSplineSurface (array) -
import OCC
import Generator as G
import Converter as C

# control points
points = G.cart((0,0,0), (1,1,1), (4,4,1))
points[1][2, 1 + 1*4] = 1.
points[1][2, 1 + 2*4] = 1.
points[1][2, 2 + 1*4] = 1.
points[1][2, 2 + 2*4] = 1.

hook = OCC.createEmptyCAD()

# from control points with uniform param
OCC._addSplineSurface(hook, points, 3)

OCC.writeCAD(hook, "out.step")
