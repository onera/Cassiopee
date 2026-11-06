# - getFaceArea (array) -
import OCC
import KCore.test as test

h = OCC.readCAD('cube.step', 'fmt_step')

# Area of all model
test.testO(OCC.getFaceArea(h), 1)

# Area of seome faces
test.testO(OCC.getFaceArea(h, [1,2]), 2)
