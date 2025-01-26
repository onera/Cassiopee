"""PolyLine mesh generator. Extension of Generator.
"""
from . import Generator as G
__version__ = G.__version__

try: range = xrange
except: pass

#=============================================================================
# PolyLineMesher pour les i-arrays ou les BAR-arrays
#=============================================================================
def polyLineMesher(polyLine, h, yplus, density):
    """Generate a multiple mesh for a polyline.
    Usage : polyLineMesher( polyLine, h, yplus, density, extension)"""
    import Geom as D
    import Converter as C
    import math
    import Transform as T

    polyLine = C.convertArray2Tetra(polyLine)
    polyLine = G.close(polyLine)

    addFactor = 0.2
    if len(polyLine) != 4:
        raise TypeError("polyLineMesher: requires a BAR array.")
    else:
        if polyLine[3] != 'BAR':
            raise TypeError("polyLineMesher: requires a BAR array.")

    f = polyLine[1]; c = polyLine[2]; ne = c.shape[1]

    # Detection de la hauteur maximum admissible
    for i in range(ne):
        ind1 = c[0,i]-1; ind2 = c[1,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        if (h > 0.9*l):
            h = 0.9*l
            print("Warning: height changed to", h,"...")
            print("...because length of line segment", i, "is", l)

    # Detection de la densite minimum
    nj = int(h*density)+1
    if (nj < 4):
        density = 4./h
        print("Warning: density changed to", density,"to have 4 points in height.")
    for i in range(ne):
        ind1 = c[0,i]-1; ind2 = c[1,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        ni = int(l*density)+1
        if (ni < 4):
            density = 4./l
            print("Warning: density changed to", density," to have 4 points in segment",i)

    # Calcul automatique de l'extension
    extension = max(int(h*density)+1, 3)

    # Calculs prealables
    curvature = D.getCurvatureAngle(polyLine)
    polyLine2 = T.addkplane(polyLine)
    np = G.getNormalMap(polyLine2)
    np = C.center2Node(np)
    np = C.normalize(np, ['sx','sy','sz'])

    nj = int(h*density)+1
    distrib = G.cart((0,0,0), (1./nj,1,1), (nj+1,1,1))
    add = max(nj * h / (20 * yplus), 1)
    add = min(int(add), 3*nj)
    add = int(addFactor*add)
    distrib = G.enforcePlusX(distrib, yplus/h, nj-1, add)
    nj = distrib[2]-1
    delta = C.array('d',nj,1,1)
    for i in range(nj):
        delta[1][0,i] = h* distrib[1][0,i]
    mesh = []; walls = []

    pool = 0; poolFirst = 0

    # Generation des maillages
    for i in range(ne):
        ind1 = c[0,i]-1; ind2 = c[1,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        ni = int(l*density)+1
        hi = extension * 1. / (ni-1)
        # ext = 1 ; extension chimere
        # ext = 0 ; TFI paroi
        # ext = -1; TFI MD + TTM
        ext1 = 1 ; ext2 = 1

        if (curvature[1][0,ind1] >= 180): # extension
            ext1 = 1
        elif (curvature[1][0,ind1] < 90+38): # TFI 2 parois
            ext1 = 0
        elif (curvature[1][0,ind1] >= 90+38): # TFI MD + TTM
            ext1 = -1
        if (curvature[1][0,ind2] >= 180): # extension
            ext2 = 1
        elif (curvature[1][0,ind2] < 90+38): # TFI 2 parois
            ext2 = 0
        elif (curvature[1][0,ind2] >= 90+38): # TFI MD + TTM
            ext2 = -1

        # Toujours extension aux frontieres
        if (ext1 != 1):
            indv = findNeighbourIndex(polyLine, ind1+1, i+1)
            if (indv == -1): ext1 = 1
        if (ext2 != 1):
            indv = findNeighbourIndex(polyLine, ind2+1, i+1)
            if (indv == -1): ext2 = 1
        #print 'mesh no=',i,' ext1=',ext1,' ext2=',ext2

        lext1 = max(0,ext1) ; lext2 = max(0,ext2)
        px1 = x1 - lext1*(x2-x1)*hi
        py1 = y1 - lext1*(y2-y1)*hi
        pz1 = z1 - lext1*(z2-z1)*hi
        px2 = x2 + lext2*(x2-x1)*hi
        py2 = y2 + lext2*(y2-y1)*hi
        pz2 = z2 + lext2*(z2-z1)*hi
        ni = ni + lext1*extension + lext2*extension
        ni = int(ni)
        line = D.line((px1,py1,pz1), (px2,py2,pz2), ni)

        if (ext1 == 1 and ext2 == 1):
            m = G.addNormalLayers(line, delta)
            m = T.reorder(m,(-1,-3,2))
        else:
            line = T.addkplane(line)
            n = G.getNormalMap(line)
            n = C.center2Node(n)
            n = C.normalize(n, ['sx','sy','sz'])

            if (ext1 == 0):
                indv = findNeighbourIndex(polyLine, ind1+1, i+1)
                indv = indv-1
                x4 = f[0,indv]; y4 = f[1,indv]; z4 = f[2,indv]
                no = math.sqrt( (x4-px1)*(x4-px1) + (y4-py1)*(y4-py1) + (z4-pz1)*(z4-pz1))
                hp = h / math.sin( curvature[1][0,ind1] * math.pi / 180. )
                px4 = px1 + hp * (x4-px1)/no
                py4 = py1 + hp * (y4-py1)/no
                pz4 = pz1 + hp * (z4-pz1)/no
            elif (ext1 == 1):
                px4 = px1 + n[1][0,0] * h
                py4 = py1 + n[1][1,0] * h
                pz4 = pz1 + n[1][2,0] * h
            elif (ext1 == -1):
                px4 = px1 + np[1][0,i] * h
                py4 = py1 + np[1][1,i] * h
                pz4 = pz1 + np[1][2,i] * h

            if (ext2 == 0):
                indv = findNeighbourIndex(polyLine, ind2+1, i+1)
                indv = indv-1
                x3 = f[0,indv]; y3 = f[1,indv]; z3 = f[2,indv]
                no = math.sqrt( (x3-px2)*(x3-px2) + (y3-py2)*(y3-py2) + (z3-pz2)*(z3-pz2))
                hp = h / math.sin( curvature[1][0,ind2] * math.pi / 180. )
                px3 = px2 + hp * (x3-px2)/no
                py3 = py2 + hp * (y3-py2)/no
                pz3 = pz2 + hp * (z3-pz2)/no
            elif (ext2 == 1):
                px3 = px2 + n[1][0,0] * h
                py3 = py2 + n[1][1,0] * h
                pz3 = pz2 + n[1][2,0] * h
            elif (ext2 == -1):
                px3 = px2 + np[1][0,i+1] * h
                py3 = py2 + np[1][1,i+1] * h
                pz3 = pz2 + np[1][2,i+1] * h

            d1 = D.line( (px1,py1,pz1), (px2,py2,pz2) )
            d2 = D.line( (px4,py4,pz4), (px3,py3,pz3) )
            d3 = D.line( (px2,py2,pz2), (px3,py3,pz3) )
            d4 = D.line( (px1,py1,pz1), (px4,py4,pz4) )
            r = G.cart((0,0,0), (1./(ni-1),1,1), (ni,1,1))
            l = D.getLength(d1)
            if (ext1 == 0):
                r = G.enforcePlusX(r, yplus/l, min(nj-1,ni-1), add)
            if (ext2 == 0):
                r = G.enforceMoinsX(r, yplus/l, min(nj-1,ni-1), add)
            r1 = G.map(d1, r)
            r2 = G.map(d2, r)
            r3 = G.map(d3, distrib)
            r4 = G.map(d4, distrib)
            r1 = T.reorder(r1, (-1,2,3))
            r2 = T.reorder(r2, (-1,2,3))
            m = G.TFI([r1, r2, r3, r4])
            #m = G.TTM(r1, r2, r3, r4, 10)
            m = T.reorder(m, (2,1,3))


        m = T.addkplane(m)
        mesh.append(m)
        # Because of reorder:
        ext = ext1; ext1 = ext2; ext2 = ext;
        # Walls
        rangesw = []
        if ( ext1 != 1 ): i1 = 1
        else: i1 = extension+1
        if ext2 != 1: i2 = m[2]
        else: i2 = m[2]-extension
        wrange = [i1,i2, 1, 1, 1, m[4]]
        rangesw.append(wrange)

        if (ext1 == 0):
            wrange = [1, 1, 1, m[3], 1, m[4]]
            rangesw.append(wrange)

        if (ext2 == 0):
            wrange = [m[2], m[2], 1, m[3], 1, m[4]]
            rangesw.append(wrange)
        walls.append(rangesw)

    return [mesh, walls, h, density]

#=============================================================================
# Return the index of neighbour of ind of element e on a BAR
def findNeighbourIndex(polyLine, index, e):
    c = polyLine[2]
    ne = c.shape[1]
    for i in range(ne):
        ind1 = c[0,i]
        ind2 = c[1,i]
        if (i != e-1):
            if (ind1 == index):
                return ind2
            if (ind2 == index):
                return ind1
    return -1
