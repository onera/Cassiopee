# - getFaceArea (array) -
import OCC

h = OCC.readCAD('cube.step', 'fmt_step')

# Area of all model
print(OCC.getFaceArea(h))

# Area of seome faces
print(OCC.getFaceArea(h, [1,2]))
