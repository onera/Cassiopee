# - splitFaces (array) -
import OCC

hook = OCC.occ.readCAD("cube.step", "fmt_step")
OCC.occ.splitFaces(hook, 20.)
OCC.occ.writeCAD(hook, "out.step", "fmt_step")
