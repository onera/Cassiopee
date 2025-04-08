# - translate (pyTre) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._translate(hook, (1,0,0))
OCC.writeCAD(hook, 'out.step', 'fmt_step')