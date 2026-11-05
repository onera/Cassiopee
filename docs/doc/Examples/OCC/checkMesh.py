# - checkMesh (array) -
import OCC.Check

hook = OCC.readCAD("cube.step", "fmt_step")
OCC.Check.checkMesh(hook)
