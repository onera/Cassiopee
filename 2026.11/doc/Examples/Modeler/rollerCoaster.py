# - Roller Coaster (array) -
import Modeler
import Converter as C
import Geom as D
import Transform as T
import Generator as G
import KCore.Vector as Vector

# COnstruit une traverse (suivant x)
b = Modeler.box((-0.5,0,0),(0.5,0.1,0.05),chamfer=0.01)

# Construit la driving line (struct 1D)
l = D.circle((-2,0,0),2.,N=30,tetas=0.,tetae=180.)
l2 = D.circle((2,0,0),2.,N=30,tetas=180.,tetae=360.)
l3 = D.line((4,0,0),(4,4,-4),N=30)
l4 = D.circle((2,4,-4),2.,N=30,tetas=0.,tetae=180.)
l5 = D.line((0,4,-4),(0,-4,-4),N=30)
l6 = D.circle((-2,-4,-4),2.,N=30,tetas=180.,tetae=360.)
l7 = D.line((-4,-4,-4),(-4,0,0),N=30)
#C.convertArrays2File([l,l2,l3,l4,l5,l6,l7],'line.plt')
l = T.join([l,l2,l3,l4,l5,l6,l7])

# Met des traverses le long de la ligne l
all = [l]
n = l[2]
coord = l[1]
e2p = None

for i in range(n-1):
    # Pt i
    Pi = [coord[0,i],coord[1,i],coord[2,i]]
    # Pt i+1
    Pip = [coord[0,i+1],coord[1,i+1],coord[2,i+1]]
    # vecteur deplacement (modele initial en 0)
    v = Pi
    # vecteur e1 (ex: y)
    e1 = Vector.sub(Pip,Pi)
    e1 = Vector.normalize(e1)
    # vecteur e2 (ex: x + beta y)
    # intersection du plan normal a e1 avec le plan x,y
    if abs(e1[1]) > 1.e-12:
        alpha = 1.
        beta = -e1[0]/e1[1]
    else:
        alpha = -e1[1]/e1[0]
        beta = 1.

    x = [1,0,0]; y = [0,1,0]
    x = Vector.mul(-alpha,x)
    y = Vector.mul(-beta,y)
    e2 = Vector.add(x,y)
    e2 = Vector.normalize(e2)
    if e2p is not None:
        if Vector.dot(e2,e2p) < -0.9:
            e2 = Vector.mul(-1,e2)
    e2p = e2
    e3 = Vector.cross(e2,e1)
    # rotate
    b2 = T.rotate(b, (0,0,0), ((1,0,0),(0,1,0),(0,0,1)), (e2,e1,e3))
    b2 = T.translate(b2, v)
    all.append(b2)

#C.convertArrays2File(all, 'out.plt')
#import sys; sys.exit()

# Construit le path ou Pts pour la camera
Pts = []
path = G.cart((0,0,0),(1,1,1),(n,1,1))

for i in range(n):
    Pts.append((coord[0,i],coord[1,i],coord[2,i]+0.2))
    path[1][0,i] = coord[0,i]
    path[1][1,i] = coord[1,i]
    path[1][2,i] = coord[2,i]+0.2

# Refine path pour l'acceleration + calcul incEye pour le eye
path = D.uniformize(path, N=2000, sharpAngle=20.)
#C.convertArrays2File(path, 'path.plt')
#import sys; sys.exit()

# Parametre les accelerations
D.setF(path, 0, 0.2) # start
D.setF(path, 400, 1.5) # mid
D.setF(path,565, 0.5) # Debut Descente
D.setF(path,818, 5.) # Fin Descente
D.setF(path,1013, 3.) # Descelaration
D.setF(path,1744, 0.4) # debut montee
D.setF(path,1999, 0.2) # fin montee
path = D.enforceh(path, h=0.0008)
#C.convertArrays2File(path, 'path.plt')

# Parametre la view incEye en fonction de la vitesse
path = C.initVars(path, 'incEye', 50.)
# incView inverse de abscisse curviligne?
s = D.getCurvilinearAbscissa(path)[1]
f = path[1]
for i in range(path[2]-1):
    path[1][3,i] = 1./(s[0,i+1]-s[0,i])*0.003
    path[1][3,i] = min(path[1][3,i],800)
    path[1][3,i] = max(path[1][3,i],20)

import CPlot
import time
CPlot.display(all, bgColor=1, displayInfo=0, mode='solid')

# MoveCamera avec ensemble de points
#for i in range(6000):
#    #time.sleep(0.001)
#    p = int(i/2000)
#    pos = i -p*2000
#    CPlot.moveCamera(Pts, N=2000, pos=pos, moveEye=True)

# MoveCamera avec array path
N = path[2]
for i in range(N*3):
    time.sleep(0.001)
    p = int(i/N)
    pos = i-p*N
    CPlot.moveCamera(path, N=N, pos=pos, moveEye=True)
    CPlot.setState(message='incEye=%d'%path[1][3,pos])
