# - reverse (array) -
import OCC

# reverse shape
hook = OCC.createEmptyCAD()
OCC._addBox(hook, (0,0,0), 1., 1., 1.)
OCC._reverse(hook)
#OCC.writeCAD(hook, "out.step")

# reverse a face
OCC._reverse(hook, [1], [], method=0)
OCC.writeCAD(hook, "out2.step")
