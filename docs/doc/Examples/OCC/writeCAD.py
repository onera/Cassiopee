# - writeCAD (array) -
import OCC

hook = OCC.occ.readCAD("cube.step", "fmt_step")
OCC.occ.writeCAD(hook, "out.step", "fmt_step")
