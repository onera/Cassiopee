# - meshTRIHO (array) -
import OCC
import Converter as C

m = OCC.meshTRIHO("cube.step", "fmt_step")
for i in m: print(i[3])
