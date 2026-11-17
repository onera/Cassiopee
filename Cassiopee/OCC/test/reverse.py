# - reverse (array) -
import OCC

# reverse shape
hook = OCC.createEmptyCAD()
OCC._addBox(hook, (0,0,0), 1., 1., 1.)
OCC._reverse(hook)
OCC.writeCAD(hook, "out.step")

# reverse a face
OCC._reverse(hook, [1], [], method=0)
OCC.writeCAD(hook, "out2.step")

# reverse an edge
hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (1,0,0))
OCC._reverse(hook, [], [1], method=0)
OCC.writeCAD(hook, 'out3.step')
