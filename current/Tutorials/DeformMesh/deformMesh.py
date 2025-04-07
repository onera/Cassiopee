# - deformation du maillage par rapport au maillage de reference -
import Converter as C
import Generator as G
import Transform as T
import Geom as D
import CPlot
import numpy

# Profil a deformer
cyl = D.circle( (0,0,0), 1. )
s = D.getCurvilinearAbscissa(cyl) 
cyl = T.addkplane(cyl) ; cyl = T.addkplane(cyl)

# Maillage a deformer
NI = cyl[2] ; NJ = 30
cylmesh = G.cylinder((0,0,0), 1., 4., 360., 0, 1., (NI,NJ,1))
distrib = G.cart((0,0,0), (1./(NJ-1),1,1), (NJ,1,1))
distrib = G.enforcePlusX(distrib, 0.001, (20,20)) 
cylmesh = G.map(cylmesh, distrib, 2)
NJ = cylmesh[3]
# Je le split pour tester le multi-domaine
cylmesh1 = T.subzone(cylmesh, (1, 1, 1), (NI//3,NJ,1))
cylmesh2 = T.subzone(cylmesh, (NI//3, 1, 1), (2*NI//3,NJ,1))
cylmesh3 = T.subzone(cylmesh, (2*NI//3, 1, 1), (NI,NJ//2,1))
cylmesh4 = T.subzone(cylmesh, (2*NI//3, NJ//2, 1), (NI,NJ,1))
cylmesh = [cylmesh1, cylmesh2, cylmesh3, cylmesh4]
cylmesh = T.addkplane(cylmesh) ; cylmesh = T.addkplane(cylmesh)

time = 0.
for i in range(200):
    # Defo de la surface (suivant une formule analytique)
    defo = C.initVars(s, '{df}=0.3*cos(%f)*sin(10*pi*{s}+%f)'%(time,time))
    defo = C.extractVars(defo, ['df'])
    defo = T.addkplane(defo) ; defo = T.addkplane(defo)
    cyl2 = T.deformNormals(cyl, defo) ; cyl2 = G.close(cyl2, 1.e-2)
    
    # Defo du maillage (suivant la surface)
    cyld = C.copy(cyl2) ; cyld[0] = 'x1,y1,z1' 
    cyl3 = C.addVars([cyl, cyld])
    cyl3 = C.initVars(cyl3, 'dx={x1}-{x}')
    cyl3 = C.initVars(cyl3, 'dy={y1}-{y}')
    cyl3 = C.initVars(cyl3, 'dz={z1}-{z}')
    cyl3 = C.extractVars(cyl3, ['x','y','z','dx','dy','dz'])
    cylmesh2 = T.deformMesh(cylmesh, cyl3, beta=7.)

    # On verifie que l'on a pas de volume negatif
    vol = G.getVolumeMap(cylmesh2)
    minvol = C.getMinValue(vol, 'vol')
    if minvol <= 0: print('BEWARE!')

    # Display
    CPlot.display([cyl2]+cylmesh2, dim=2)
    time += 0.05

C.convertArrays2File([cyl2]+cylmesh2, 'out.plt')
