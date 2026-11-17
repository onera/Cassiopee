# - splitFace (array) -
import OCC

# split planar face
hook = OCC.createEmptyCAD()
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)
OCC._fillHole(hook, [1], 0)
OCC._splitFace(hook, 1, (25,25,0), (-25,-25,0))
#OCC._splitFace(hook, 1, (0.5,0.5,0), (-0.5,-0.5,0))
OCC.writeCAD(hook, "out.step")

# split non planar face
hook = OCC.createEmptyCAD()
OCC._addSphere(hook, (0,0,0), 1.)
OCC._splitFace(hook, 1, (-50,-50,0), (50,50,0))
OCC.writeCAD(hook, "out2.step")

# split non planar face
hook = OCC.createEmptyCAD()
OCC._addLine(hook, (0,0,0), (10,0,0))
OCC._addLine(hook, (10,0,0), (10,10,3))
OCC._addLine(hook, (10,10,3), (0,10,0))
OCC._addLine(hook, (0,10,0), (0,0,0))
OCC._fillHole(hook, [1,2,3,4], [], 1)
OCC._splitFace(hook, 1, (-50,-50,0), (50,50,0))
OCC.writeCAD(hook, "out3.step")
