# - splitFace (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._splitFace(hook, 1, (25,25,-50), (-25,25,0))
OCC.writeCAD(hook, "out.step", "fmt_step")
