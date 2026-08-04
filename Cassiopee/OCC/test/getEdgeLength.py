# - getEdgeLength (array) -
import OCC

h = OCC.readCAD('cube.step', 'fmt_step')

# Length of all edges
print(OCC.getEdgeLength(h))

# Length of one edge
print(OCC.getEdgeLength(h, [1]))
