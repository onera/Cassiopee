# - removeFaces (array) -
import OCC

hook = OCC.readCAD("cube.step")

OCC._extractFaces(hook, [1,2])

OCC.writeCAD(hook, 'out.step')
