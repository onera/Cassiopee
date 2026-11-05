# - getFaceVolume (array) -
import OCC
import KCore.test as test

h = OCC.readCAD('cube.step', 'fmt_step')

# Volume of all model
vol = OCC.getFaceVolume(h)
test.testO(vol, 1)

# Volume of some faces (must be closed)
vol = OCC.getFaceVolume(h, [1,2,3,4,5,6])
test.testO(vol, 2)
