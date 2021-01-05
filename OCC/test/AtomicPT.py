# - Atomic functions test -
import OCC
import Converter as C
hook = OCC.occ.readCAD("hammer.iges", "fmt_iges")
out = OCC.occ.meshGlobalEdges(hook, 10.)

C.convertArrays2File(out, 'out.plt')