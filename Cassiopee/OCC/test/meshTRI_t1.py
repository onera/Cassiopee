# - meshTRI (array) -
import OCC
import KCore.test as test

m = OCC.meshTRI("cube.step", "fmt_step")
test.testA(m, 1)