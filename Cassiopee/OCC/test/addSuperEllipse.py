# - addSuperEllipse (array) -
import OCC

hook = OCC.createEmptyCAD()
# add squircle
OCC._addSuperEllipse(hook, (0,0,0), 1., 1., 4, 36, False)

# add super ellipse
OCC._addSuperEllipse(hook, (0,3,0), 1., 1.5, 4, 36, False)

OCC.writeCAD(hook, "out.step")
