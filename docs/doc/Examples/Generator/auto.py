#==============================================================================
# - The *Cassiopee* Project -
# Template file for automatic mesh generation
# (c) 2005-2009 by ONERA
# Input : a geometry 1D curve
# Output : 2D meshes
# The 1D curve is split following curvature
# Then hyperbolic meshes are generated
#==============================================================================

#==============================================================================
# User data
#==============================================================================
# Geometry file
GEOM_FILE = "lines.tp"
# Total height of mesh
MSH_HEIGHT = 0.1
# First cell height
BND_HEIGHT = 0.001
# Density of points (number of points for a length of 1) in i and j directions
DI = 30
DJ = 200

#==============================================================================
import Transform as T
import Converter as C
import Generator as G
import Geom as D

a = C.convertFile2Arrays('lines.tp', 'fmt_tp')
geom = a[0]

# Split geom curve
mshList = T.split( geom, 1.e-8 )
n = len( mshList )

c = 0
for x in mshList:

    # Define the 1D distributions
    l = D.getLength( x )
    ni = int(DI*l)
    nj = int(DJ*MSH_HEIGHT)
    hi = 1./(ni-1)
    hj = MSH_HEIGHT/(nj-1)
    hk = 1.
    distrib = G.cart((0.,0.,0.), (hi, hj, hk), (ni, nj, 1))
    distrib = G.enforcePlusY(distrib, BND_HEIGHT, (10,10))
    # Recompute distrib data
    ni = distrib[2]
    nj = distrib[3]
#   if (hj > hi*l):
#       print("Warning, mesh generator may be unstable because")
#       print("cell height is greater than cell length")
    temp1 = x
    angleLeft = 180.
    angleRight = 180.

    if (c > 0):
        # Left extension
        ln1 = distrib[1][0][1]; ln2 = distrib[1][0][0];
        dl = ln1-ln2
        dl = dl * l
        distrib = G.addPointInDistribution(distrib, 1)
        niv = mshList[c-1][2]
        temp2 = mshList[c-1]
        dn = D.getDistantIndex(temp2, niv, -dl)
#        print('dn gauche ', dn, dl, niv)
        temp2 = T.subzone( temp2, (dn,1,1), (niv,1,1) )
        temp1 = T.join( temp2, temp1 )
        curv = D.getCurvatureAngle( temp1 )
        angleLeft = curv[1][0][niv-dn+1]

    if (c < n-1):
        # Right extension
        ln1 = distrib[1][0][ni-1]; ln2 = distrib[1][0][ni-2];
        dl = ln1-ln2
        dl = dl * l
        distrib = G.addPointInDistribution(distrib, ni)

        temp2 = mshList[c+1]
        niv = temp2[2]
        dn = D.getDistantIndex(temp2, 1, dl)
#        print('dn droit ', dn, dl)
        temp2 = T.subzone( temp2, (1,1,1), (dn,1,1) )
        temp2 = T.join( temp1, temp2 )
        temp1 = temp2
        curv = D.getCurvatureAngle( temp1 )
        angleRight = curv[1][0][temp1[2]-dn+1]

    print("Angles =", angleLeft, " ", angleRight)
    temp1 = G.hyper2D3(temp1, distrib, "C", angleLeft, angleRight)
    C.convertArrays2File([temp1], "mesh"+repr(c)+".plt", "bin_tp")
    c += 1
