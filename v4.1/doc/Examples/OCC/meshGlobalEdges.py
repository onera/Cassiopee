# - meshGlobalEdges (array) -
import OCC
import Converter as C

hook = OCC.occ.readCAD("cube.step", "fmt_step")
edges = OCC.occ.meshGlobalEdges1(hook, 10.)
C.convertArrays2File(edges, 'out.plt')
