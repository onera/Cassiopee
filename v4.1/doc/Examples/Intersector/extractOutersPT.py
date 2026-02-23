# - extractOuterLayers (pyTree) -
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Generator.PyTree as G
t = C.convertFile2PyTree('boolNG_M1.tp')
t = C.conformizeNGon(t); t = G.close(t)
t = XOR.extractOuterLayers(t, 1, discard_external=0, output_remaining=True)

C.convertPyTree2File(t, "out.cgns")
