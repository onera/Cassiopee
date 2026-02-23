# - recoverGlobalIndex (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)

b = T.splitNParts(a, 2)
C._initVars(b[0], 'f=1')
C._initVars(b[1], 'f=2')

C._recoverGlobalIndex(a, b[0])
C._recoverGlobalIndex(a, b[1])
C._rmVars(a, 'globalIndex')

C.convertPyTree2File(a, 'out.cgns')
