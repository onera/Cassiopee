# - mergeFaces (array) -
import OCC

hook = OCC.readCAD("cube.step")
OCC._splitFaces(hook, 200.)
# merge given faces
OCC._mergeFaces(hook, [1,2,3])
# merge all faces
OCC._mergeFaces(hook)
OCC.writeCAD(hook, "out.step")
