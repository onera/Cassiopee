# - convertArrays2ZoneNode
import Converter.PyTree as C
import Generator as G

a = G.cart((0.,0.,0.),(1.,1.,1.),(5,5,5))
zone = C.convertArrays2ZoneNode('cart',[a])
print(zone)
