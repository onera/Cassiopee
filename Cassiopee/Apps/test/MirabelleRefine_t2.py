# - Mesh.Mirabelle3 -
# propagate a new edge distribution imposing the edge values
import Converter.PyTree as C
import Transform.PyTree as T
import Apps.Mesh.Mirabelle3 as Mirabelle
import KCore.test as test

t = C.convertFile2PyTree("premirabelle.cgns")

zonesToRefine = ["zone.10", "zone.9"]
dirs = [1,1]
Mirabelle._refine(t, zonesToRefine, dirs, refined={}, factor=2)
test.testT(t, 1)
