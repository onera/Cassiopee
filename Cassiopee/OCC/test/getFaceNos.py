# - getFaceNos (array) -
import OCC

hook = OCC.readCAD("cube.step", "fmt_step")

ret = OCC.getFaceNos(hook, "Part1")
print(ret)
