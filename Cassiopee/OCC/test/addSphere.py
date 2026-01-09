# - addSphere (array) -
import OCC

hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0.,0.,0.), 1.)
#OCC._addSphere(hook, (2.,0.,0.), 0.5)
OCC.writeCAD(hook, "out.step")
