# - meshGlobalEdges (array) -
import OCC
import KCore.test as test

hook = OCC.readCAD("cube.step", "fmt_step")
edges = OCC.occ.meshGlobalEdges1(hook, 10.)
test.testA(edges, 1)
