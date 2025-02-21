# - writeCAD (pyTree) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC.writeCAD(hook, "out.step", "fmt_step")
