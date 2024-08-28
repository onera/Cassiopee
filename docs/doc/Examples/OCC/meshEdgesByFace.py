# - meshEdgesByFace (array) -
import OCC
import Converter as C

hook = OCC.occ.readCAD("cube.step", "fmt_step")
edges = OCC.occ.meshEdgesByFace(hook, 1, 10, -1.)

C.convertArrays2File(edges, 'out.plt')
