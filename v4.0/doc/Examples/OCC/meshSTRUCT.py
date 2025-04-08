# - meshSTRUCT (array) -
import OCC
import Converter as C

m = OCC.meshSTRUCT("cube.step", "fmt_step")
C.convertArrays2File(m, 'out.plt')
