# - rotate (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._rotate(hook, (0,0,0), (0,0,1), 30.)
OCC._rotate(hook, (0,0,0), (0,0,1), 30., faceList=[1])
OCC.writeCAD(hook, 'out.step', 'fmt_step')