# - fixShape (array) -
import OCC

hook = OCC.readCAD("cube.step")

# fix all shape
OCC._fixShape(hook, fixShape=0, fixWires=0, unifyEdges=1, unifyFaces=1,
              tol=1.e-6, angularDeflection=1., linearDeflection=0.01)
OCC.writeCAD(hook, 'out.step')
