# - rmBCDataVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'jmax')
C._initBCDataVars(a,'{var1}=1.')
C._initBCDataVars(a,'{var2}=2.')
C._initBCDataVars(a,'{var3}=3.')
a = C.rmBCDataVars(a,'var1')
a = C.rmBCDataVars(a,['var2','var3'])

