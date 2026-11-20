# - boolean (array) -
import OCC

# union of boxes
hook = OCC.createEmptyCAD()
OCC._addBox2(hook, (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1))
OCC._addBox2(hook, (0.8,0.8,0.8), (1.2,0.8,0.8), (1.2,1.2,0.8), (0.8,1.2,0.8),
             (0.8,0.8,1.2), (1.2,0.8,1.2), (1.2,1.2,1.2), (0.8,1.2, 1.2))
OCC._boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 0, 0, 0)
OCC.writeCAD(hook, "box1.step")

# union of lofts
hook = OCC.createEmptyCAD()
OCC._addEllipse(hook, (0,0,0), (0,0,1), 1., 2.)
OCC._addEllipse(hook, (0,0,5), (0,0,1), 1., 2.)
OCC._loft(hook, [1,2], [])
OCC._fillHole(hook, [1])
OCC._fillHole(hook, [2])
hook2 = OCC.createEmptyCAD()
OCC._addEllipse(hook2, (0,0,0), (0,0,1), 1., 2.)
OCC._addEllipse(hook2, (0,0,5), (0,0,1), 1., 2.)
OCC._loft(hook2, [1,2], [])
OCC._fillHole(hook2, [1])
OCC._fillHole(hook2, [2])
OCC._rotate(hook2, (0,0,0), (1,0,0), 90.)
hook = OCC.mergeCAD([hook, hook2])
OCC._boolean(hook, [1,2,3], [4,5,6], 0, 0, 0)
OCC.writeCAD(hook, 'loft1.step')

# difference of spheres
hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0,0,0), 1.)
OCC._addSphere(hook, (1,0,0), 0.2)
#OCC._boolean(hook, [1], [2], 1, 0, 0) # FAIL
OCC.writeCAD(hook, "sphere1.step")

# difference of boxes
hook = OCC.createEmptyCAD()
OCC._addBox2(hook, (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1))
OCC._addBox2(hook, (0.8,0.8,0.8), (1.2,0.8,0.8), (1.2,1.2,0.8), (0.8,1.2,0.8),
             (0.8,0.8,1.2), (1.2,0.8,1.2), (1.2,1.2,1.2), (0.8,1.2, 1.2))
OCC._boolean(hook, [i for i in range(1,7)], [i for i in range(7,13)], 1, 0, 0)
OCC.writeCAD(hook, "box2.step")
