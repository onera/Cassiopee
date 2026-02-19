# - trimFaces (pyTree) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cubeNoTrim.step", "fmt_step")

ret = OCC.getFaceNameInOCAF(hook)
cube1 = ret[3]
cube2 = ret[5]
#OCC._translate(hook, (+300.,0,0), cube2)
OCC._trimFaces(hook, cube1, cube2)
OCC._sewing(hook, tol=1.e-2)

OCC.writeCAD(hook, 'out.step', 'fmt_step')
