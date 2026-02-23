# - scale (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._scale(hook, 0.5, (0,0,0))
OCC.writeCAD(hook, 'out.step', 'fmt_step')