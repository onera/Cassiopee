# - recoverGlobalIndex (array) -
import Converter as C
import Generator as G
import Transform as T

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)

b = T.splitNParts(a, 2)
C._initVars(b[0], 'f=1')
C._initVars(b[1], 'f=2')

C._recoverGlobalIndex(a, b[0])
C._recoverGlobalIndex(a, b[1])
a = C.rmVars(a, 'globalIndex')

C.convertArrays2File(a, 'out.plt')
