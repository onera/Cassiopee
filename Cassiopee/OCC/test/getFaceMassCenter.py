# - getFaceMassCenter (array) -
import OCC

h = OCC.readCAD('cube.step', 'fmt_step')

# Center of mass of all model
print(OCC.getFaceMassCenter(h))

# Center of mass of some faces (must be closed)
print(OCC.getFaceMassCenter(h, [1,2,3,4,5,6]))
