# - addFillet (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._addFillet(hook, [1,2,3,4,5,6,7,8,9,10,11,12], radius=10.)
OCC.writeCAD(hook, 'out.step', 'fmt_step')
