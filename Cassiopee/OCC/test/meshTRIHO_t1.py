# - meshTRIHO (array) -
import OCC
import Converter as C
import KCore.test as test

m = OCC.meshTRIHO("cube.step", "fmt_step")
for i in m: print(i[3])
test.testA(m, 1)