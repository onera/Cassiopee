# - getBoundingBox (array) -
import OCC

# object
hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0,0,0), 1.)

bbox = OCC.getBoundingBox(hook)
print(bbox)
