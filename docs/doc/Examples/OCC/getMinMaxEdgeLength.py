# - getMinMaxEdgeLength (array) -
import OCC

h = OCC.readCAD('cube.step', 'fmt_step')

# Min/Max of all edges
print(OCC.getMinMaxEdgeLength(h))

# Area of some faces
print(OCC.getMinMaxEdgeLength(h, [1,2]))
