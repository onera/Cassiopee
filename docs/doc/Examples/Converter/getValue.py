# - getValue (array) -
import Converter as C
import Generator as G

# Structured array
Ni = 40; Nj = 50; Nk = 20
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
# Get variable values contained in a (x,y,z) in point (10,1,1)
print(C.getValue(a, (10,1,1)))
#>> [0.23076923076923075, 0.0, 0.0]
print(C.getValue(a, 9)) # It is the same point!
#>> [0.23076923076923075, 0.0, 0.0]

# Unstructured array
Ni = 40; Nj = 50; Nk = 20
a = G.cartTetra((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
print(C.getValue(a, 9))
#>> [0.23076923076923075, 0.0, 0.0]
