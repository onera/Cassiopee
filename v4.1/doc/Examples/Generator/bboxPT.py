# - bbox (pyTree) -
import Generator.PyTree as G
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
print(G.bbox(a))
