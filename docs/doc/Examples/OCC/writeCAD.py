# - writeCAD (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC.writeCAD(hook, "out.step", "fmt_step")
