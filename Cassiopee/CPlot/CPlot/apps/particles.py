# - particles -
# Ajoute des particules se deplacant suivant la vitesse definie dans un champ
import Converter as C
import Post as P
import CPlot

# Solution en noeuds et en centres
a = C.convertFile2Arrays('outputn.plt')
sol = C.convertFile2Arrays('output.plt')
CPlot.display(sol, dim=2)

# Entree du point d'emission
print("Cliquez pour definir le point d'emission...")

Point = []
while Point == []:
    Point = CPlot.getActivePoint()

# Point d'emission
print("Point d'emission :", Point)

# Particules
np = 20; nq = 20
nodes = C.array('x,y,z,ro,rou,rov,row,roE,cellN', np*nq, 0, 'NODE')
x = nodes[1]
for j in range(nq):
    for i in range(np):
        x[0,i+j*np] = Point[0]+0.04*i
        x[1,i+j*np] = Point[1]+0.04*j
        x[2,i+j*np] = Point[2]

# Advection
t = 0.; dt = 0.1; tfinal = 100.

while t < tfinal:
    nodes = P.extractMesh(a, nodes)
    x = nodes[1]
    for p in range(np*nq):
        px = x[0,p] + x[4,p] * dt;
        py = x[1,p] + x[5,p] * dt;
        pz = x[2,p] + x[6,p] * dt;
        x[0,p] = px; x[1,p] = py; x[2,p] = pz;
    CPlot.display(sol+[nodes])
