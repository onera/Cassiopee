# - splitFaces (pyTree) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._splitFaces(hook, 20.)
OCC.writeCAD(hook, "out.step", "fmt_step")
