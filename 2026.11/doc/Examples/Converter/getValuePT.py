# - getValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

# Structured array
Ni = 40; Nj = 50; Nk = 20
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
# Get variable values contained in a in point (10,1,1)
print(C.getValue( a, 'CoordinateX', (10,1,1) ))
#>> 0.230769230769
print(C.getValue( a, 'CoordinateX', 9 )) # It's the same point
#>> 0.230769230769
print(C.getValue( a, 'nodes:CoordinateX', 9 )) # It's the same point
#>> 0.230769230769
print(C.getValue( a, 'GridCoordinates', 9 )) # return [x,y,z]
#>> [0.23076923076923075, 0.0, 0.0]
print(C.getValue( a, ['CoordinateX', 'CoordinateY'], 9 )) # return [x,y]
#>> [0.23076923076923075, 0.0]
# Unstructured array
Ni = 40; Nj = 50; Nk = 20
a = G.cartTetra((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
print(C.getValue( a, 'CoordinateX', 9 ))
#>> 0.230769230769
