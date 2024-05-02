# - setValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

# Structured array
Ni = 40; Nj = 50; Nk = 20
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
C.setValue(a, 'CoordinateX', (10,1,1), 0.25)
C.setValue(a, 'GridCoordinates', (11,1,1), [0.3,0.2,0.1]); print(a)

# Unstructured array
Ni = 40; Nj = 50; Nk = 20
a = G.cartTetra((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
C.setValue(a, 'CoordinateX', 9, 0.1 ); print(a)
