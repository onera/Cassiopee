# - meshGlobalEdges (array) -
import OCC
import Converter as C
import KCore.test as test

hook = OCC.occ.readCAD("cube.step", "fmt_step")
edges = OCC.occ.meshGlobalEdges1(hook, 10.)
test.testA(edges, 1)
