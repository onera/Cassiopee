# - removeFaces (array) -
import OCC

hook = OCC.readCAD("cube.step")

#faceList = OCC.getFaceNos(hook, 'Cube')
OCC._removeFaces(hook, [1,2])

OCC.writeCAD(hook, 'out.step')
