# - removeFaces (array) -
import OCC

hook = OCC.readCAD("cubeNoTrim.step", "fmt_step")

ret = OCC.getFaceNameInOCAF(hook)
cube1 = ret[3]
OCC._removeFaces(hook, cube1)

OCC.writeCAD(hook, 'out.step', 'fmt_step')
