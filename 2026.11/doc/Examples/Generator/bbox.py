# - bbox (array) -
import Generator as G
a = G.cart((0.,0.,0.),(0.1,0.1,1.),(20,20,20))
b = G.cart((12.,0.,0.),(0.1,0.1,1.),(20,20,20))
print(G.bbox(a))
print(G.bbox([a, b]))
