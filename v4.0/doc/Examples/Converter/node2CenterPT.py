# - node2Center (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (30,40,1))
a = C.initVars(a, '{Density}=2*{CoordinateX}+{CoordinateY}')

# node2Center: passe une variable en centres (dans la meme zone)
a = C.node2Center(a, 'Density')
C.convertPyTree2File(a, 'out1.cgns')

# node2Center: cree une nouvelle zone contenant les centres
a = G.cart((0,0,0), (1,1,1), (30,40,1))
a = C.initVars(a, '{Density}=2*{CoordinateX}+{CoordinateY}')
b = C.node2Center(a); b[0] = a[0]+'_centers'
C.convertPyTree2File([a,b], 'out2.cgns')
