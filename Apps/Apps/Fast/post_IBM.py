#--------------------------------------------------
# Ce script permet de calculer l'epaisseur de la couche limite pour chaque point IB.
# Si 2D : chargement du restart/tc_restart/case
# Si 3D : chargement d'une slice du restart/tc_restart/case
#--------------------------------------------------
import Converter.PyTree as C
import Post.PyTree as P
import Connector.ToolboxIBM as IBM
import Converter.Internal as Internal
import Converter.Filter as Filter
import Generator.PyTree as G
import Geom.PyTree as D
from math import *
import numpy

modulo    = 1 # Choisi de ne traiter qu'un point IB tous les modulo points
distance  = 0.1 # Longueur maximale des lignes sondes
tolerance = 1e-8 # Valeur a partir de laquelle on considere nu == 0
coord     = 2.1084 # coordY de la slice pour les cas 3D

#-----------------------------------------------------
# h = Filter.Handle('restart_slice.cgns')
h = Filter.Handle('restart.cgns')
t = h.loadSkeleton()
h._loadZonesWoVars(t)
h._loadVariables(t, var=['centers:TurbulentSANuTilde'])

tb = C.convertFile2PyTree("case.cgns")
tb = C.convertArray2Tetra(tb)

tc = C.convertFile2PyTree('tc_restart.cgns')

model = Internal.getNodeFromName(tb, 'GoverningEquations')
model = Internal.getValue(model)

snear = []
sn = Internal.getNodesFromName(tb, 'snear')
for node in sn:
    snear.append(Internal.getValue(node))
snear = numpy.min(snear)

dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
dimPb = Internal.getValue(dimPb)

if dimPb == 3:
    zw = IBM.extractIBMWallFields(tc)
    zw = P.selectCells(zw, '({CoordinateY}>%f-0.001)&({CoordinateY}<%f+0.001)'%(coord,coord))
    coordY = Internal.getNodeFromName(zw, 'CoordinateY')
    nbPts  = coordY[1].shape[0]
    Internal.setValue(coordY, numpy.ones(nbPts)*coord)
    Internal.setValue(Internal.getNodeFromName(zw, 'CoordinateY_PI'), numpy.ones(nbPts)*coord)
    zw = G.close(zw)
else:
    zw = IBM.extractIBMWallFields(tc)
    zw = G.close(zw)

#-----------------------------------------------------
C._initVars(zw,'{boundaryLayerThickness}=1000')

coordX_PW   = Internal.getNodeFromName(zw,'CoordinateX_PW')[1]
coordY_PW   = Internal.getNodeFromName(zw,'CoordinateY_PW')[1]
coordZ_PW   = Internal.getNodeFromName(zw,'CoordinateZ_PW')[1]

coordX_PI   = Internal.getNodeFromName(zw,'CoordinateX_PI')[1]
coordY_PI   = Internal.getNodeFromName(zw,'CoordinateY_PI')[1]
coordZ_PI   = Internal.getNodeFromName(zw,'CoordinateZ_PI')[1]

blThickness = Internal.getNodeFromName(zw,'boundaryLayerThickness')[1]

nIBC    = coordX_PW.shape[0]
nbPts   = int(distance/snear) #Le nombre de points de chaque ligne est calcule a partir du plus petit snear

for i in range(0,nIBC,modulo):
    print('****', i, '/',nIBC)

    # Les normes sont recalculees a partir des points wall et images
    nx = coordX_PI[i]-coordX_PW[i]
    ny = coordY_PI[i]-coordY_PW[i]
    nz = coordZ_PI[i]-coordZ_PW[i]
    if dimPb==2: nz=0
    norm = sqrt(nx*nx+ny*ny+nz*nz)

    PW = (coordX_PW[i],coordY_PW[i],coordZ_PW[i])
    # 0.0005 est present pour se decaler legerement de la paroi
    # Ce nouveau point partant de PW est appele P1
    P1 = (PW[0]+nx/norm*0.0005, PW[1]+ny/norm*0.0005, PW[2]+nz/norm*0.0005) # attention au sens de la normale
    P2 = (P1[0]+nx/norm*distance, P1[1]+ny/norm*distance, P1[2]+nz/norm*distance)

    if i==0: line = D.line(P1, P2, nbPts); line = C.convertArray2Hexa(line) #Initialisation de la zone line
    else:
        # On ajoute chaque nouvelle linelet a line
        # On ajoute a la main les coordonnees et elements de chaque linelet dans la zone globale
        linelet = D.line(P1, P2, nbPts); linelet = C.convertArray2Hexa(linelet)
        coordX = Internal.getNodeFromName(line, 'CoordinateX')
        coordY = Internal.getNodeFromName(line, 'CoordinateY')
        coordZ = Internal.getNodeFromName(line, 'CoordinateZ')

        eltRa = Internal.getNodeFromName(line, 'ElementRange')
        eltCo = Internal.getNodeFromName(line, 'ElementConnectivity')

        Internal.setValue(coordX, numpy.concatenate((coordX[1], Internal.getNodeFromName(linelet, 'CoordinateX')[1]), axis=0))
        Internal.setValue(coordY, numpy.concatenate((coordY[1], Internal.getNodeFromName(linelet, 'CoordinateY')[1]), axis=0))
        Internal.setValue(coordZ, numpy.concatenate((coordZ[1], Internal.getNodeFromName(linelet, 'CoordinateZ')[1]), axis=0))

        Internal.setValue(eltRa, numpy.array((1,eltRa[1][1]+Internal.getNodeFromName(linelet, 'ElementRange')[1][1]), dtype='int32'))
        Internal.setValue(eltCo, numpy.concatenate((eltCo[1], Internal.getNodeFromName(linelet, 'ElementConnectivity')[1]+eltCo[1][-1]), axis=0))

        # On n'oublie pas de mettre a jour les dimensions de line
        line[1][0,0] = eltCo[1][-1]
        line[1][0,1] = eltRa[1][1]

coordX = Internal.getNodeFromName(line, 'CoordinateX')[1]
coordY = Internal.getNodeFromName(line, 'CoordinateY')[1]
coordZ = Internal.getNodeFromName(line, 'CoordinateZ')[1]

# Unique extractmesh a partir du restart
line = P.extractMesh(t, line)

C.convertPyTree2File(line, "line.cgns")

# On boucle sur chaque ligne pour trouver le point ou s'annule nu
# On ajoute ensuite la distance de couche limite obtenue dans le fichier wall final
nu = Internal.getNodeFromName(line, 'TurbulentSANuTilde')[1]
for i in range(0,nIBC,modulo):
    for j in range(i/modulo*nbPts+1, i/modulo*nbPts+nbPts):
        if nu[j] < tolerance:
            blThickness[i] = sqrt((coordX[j]-coordX_PW[i])**2 + (coordY[j]-coordY_PW[i])**2 + (coordZ[j]-coordZ_PW[i])**2 )
            break

C.convertPyTree2File(zw, 'wall.cgns')
