# - mergeFaces (array) -
import OCC.PyTree as OCC

hook = OCC.readCAD("cube.step", "fmt_step")
OCC._splitFaces(hook, 200.)
# merge given faces
OCC._mergeFaces(hook, [1,2,3])
# merge all faces
OCC._mergeFaces(hook)
OCC.writeCAD(hook, "out.step", "fmt_step")
