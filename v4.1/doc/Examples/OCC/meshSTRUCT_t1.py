# - meshSTRUCT (array) -
import OCC
import Converter as C
import KCore.test as test

m = OCC.meshSTRUCT("cube.step", "fmt_step")
test.testA(m, 1)