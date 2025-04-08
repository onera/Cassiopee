# - exportAccumulatorPerPsi (pyTree) -
import Post.Rotor as PR
import Converter.PyTree as C

accu = {}
for psi in range(0, 360, 10):
    for rad in range(0, 10):
        accu[(psi,rad)] = [psi+rad*0.1]

z = PR.exportAccumulatorPerPsi(accu, psi=10., vars=['F'])
C.convertPyTree2File(z, 'out.cgns')