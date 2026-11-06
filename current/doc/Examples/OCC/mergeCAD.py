# - mergeCAD (array) -
import OCC

hook1 = OCC.readCAD('cube.step')

hook2 = OCC.createEmptyCAD()
OCC._addSphere(hook2, (0,0,0), 100., name='sphere1')

#OCC._mergeCAD([hook1, hook2])
#OCC.writeCAD(hook1, "out.step")

hook = OCC.mergeCAD([hook1, hook2])
OCC.writeCAD(hook, 'out.step')
