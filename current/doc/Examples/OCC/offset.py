# - offset (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._offset(hook, 10.)
OCC.writeCAD(hook, 'out.step', 'fmt_step')
