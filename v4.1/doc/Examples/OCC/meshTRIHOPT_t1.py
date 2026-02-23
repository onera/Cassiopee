# - meshTRIHO (pyTree) -
import OCC.PyTree as OCC
import Converter.PyTree as C
import KCore.test as test

m = OCC.meshTRIHO("cube.step", "fmt_step")
test.testT(m, 1)