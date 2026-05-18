# - getFaceMassCenter (array) -
import OCC
import KCore.test as test

h = OCC.readCAD('cube.step', 'fmt_step')

# Center of mass of all model
Xc = OCC.getFaceMassCenter(h)
test.testO(Xc, 1)

# Center of mass of some faces (must be closed)
Xc = OCC.getFaceMassCenter(h, [1,2,3,4,5,6])
test.testO(Xc, 2)
