# - subzone (array) -
import Converter as C
import Transform as T
import Generator as G

# Structure
a = G.cart((0,0,0), (1,1,1), (10,20,10))
a = T.subzone(a, (3,3,3), (7,8,5))

# Structure avec indices negatif
e = G.cart((0,0,0), (1,1,1), (10,20,10))
e = T.subzone(e, (1,1,1), (-1,-1,-2)) # imax,jmax,kmax-1

# Non structure. Indices de noeuds -> retourne elts
b = G.cartTetra((0,0,0), (1,1,1), (2,2,2))
b = T.subzone(b, [2,6,8,5])

# Non structure. Indices d'elements -> Retourne elts
# Les indices d'elements commencent a 0
c = G.cartTetra((0,0,0), (1,1,1), (5,5,5))
c = T.subzone(c, [0,1], type='elements')

# Non structure. Indices de faces:
# Pour les maillages BAR, TRI, TETRA... indFace=indElt*nbreFaces+noFace
# les noFace commence a 1
# Pour les NGONS... indices des faces
# -> Retourne les faces
d = G.cartTetra((0,0,0), (1,1,1), (2,2,2))
d = T.subzone(d, [1,2,3], type='faces')

C.convertArrays2File([a,b,c,d,e], 'out.plt')
