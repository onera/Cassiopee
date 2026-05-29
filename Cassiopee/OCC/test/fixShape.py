# - fixShape (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._fixShape(hook, fixShape=1, fixWires=1, unify=1, tol=1.e-6)
OCC.writeCAD(hook, 'out.step', 'fmt_step')
