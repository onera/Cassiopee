# - splitMultiplePts (array) -
import Generator as G
import Transform as T
import Converter as C

z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,1))
z1 = T.subzone(z0,(1,1,1),(5,10,1))
z2 = T.subzone(z0,(5,1,1),(10,5,1))
z3 = T.subzone(z0,(5,5,1),(10,10,1))
zones = [z1,z2,z3]
zones = T.splitMultiplePts(zones,dim=2)
C.convertArrays2File(zones, 'out.plt')
