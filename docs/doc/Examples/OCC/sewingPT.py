# - sewing (pyTree) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._sewing(hook, tol=1.e-6, listFaces=[1,2])
OCC._sewing(hook, tol=1.)
OCC.writeCAD(hook, 'out.step', 'fmt_step')
