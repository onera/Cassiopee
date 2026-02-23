# - rotate (pyTree) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._rotate(hook, (0,0,0), (0,0,1), 30.)
OCC.writeCAD(hook, 'out.step', 'fmt_step')