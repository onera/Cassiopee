# - getFaceVolume (array) -
import OCC

h = OCC.readCAD('cube.step', 'fmt_step')

# Volume of all model
print(OCC.getFaceVolume(h))

# Volume of some faces (must be closed)
print(OCC.getFaceVolume(h, [1,2,3,4,5,6]))
