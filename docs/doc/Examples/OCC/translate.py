# - translate (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._translate(hook, (1,0,0))
OCC._translate(hook, (0,5,0), listFaces=[1])
OCC.writeCAD(hook, 'out.step', 'fmt_step')