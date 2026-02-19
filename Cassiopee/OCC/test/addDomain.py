# - addDomain (array) -
import OCC

# add box around sphere
hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0,0,0), 1.)
OCC._addDomain(hook, type="box", dfar=(10,15,10))
OCC.writeCAD(hook, 'out.step')

# add half sphere around sym box
hook = OCC.createEmptyCAD()
OCC._addBox(hook, (0,0,0), 0.5, 1., 1.)
OCC._removeFaces(hook, faceList=[3])
OCC._addDomain(hook, type="half-sphere", dfar=10, plane="xmin")
OCC.writeCAD(hook, 'out.step')
