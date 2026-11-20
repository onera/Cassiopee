# - addEllipse (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addEllipse(hook, (0,0,0), (0,0,1), 2., 1.)
OCC._addEllipse(hook, (3,0,0), (0,0,1), 1., 2.)
OCC.writeCAD(hook, "out.step")
