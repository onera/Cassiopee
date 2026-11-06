# - exportAccumulatorMap (pyTree) -
import Post.Rotor as PR
import Converter.PyTree as C

accu = {}
for psi in range(0, 365, 10):
    for rad in range(0, 10):
        accu[(psi,rad)] = [psi+rad*0.1]

z = PR.exportAccumulatorMap(accu, vars=['F'])
C.convertPyTree2File(z, 'out.cgns')
