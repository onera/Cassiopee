"""PolyC1 mesh generator. Extension of Generator.
"""
from . import Generator as G
import math
try:
    import Geom as D
    import Post as P
    import Converter as C
    import Transform as T
    import numpy as N
except: raise ImportError("PolyC1 module requires numpy, Geom, Post, Transform and Converter modules.")
__version__ = G.__version__

from Converter.Internal import E_NpyInt

try: range = xrange
except: pass

#=============================================================================
# PolyC1Mesher pour les i-arrays
#=============================================================================
def polyC1Mesher(curve, h, yplus, density, splitCrit=10., dalpha=5.,depth=1):
    """Generate a multiple mesh for a  C1 curve (polyC1).
    Usage: polyC1Mesher(curve, h, yplus, density, splitCrit, dalpha, depth)"""
    # Determination de l'angle a partir duquel on construit un raccord
    # coincident au lieu d'une extension recouvrante
    # Evite de plus les extrapolations locales pour les angles tres obtus
    alpha0 = 4.*math.atan((2*depth-1)*yplus*density)*180/math.pi
    alpha0 = max(alpha0, dalpha)
    alphaMax = 360.-alpha0

    # Split curve
    curve = G.close(curve)    
    curves = T.splitCurvatureRadius(curve, splitCrit)
    #C.convertArrays2File(curves, 'split.plt')

    addFactor = 0.2

    eps = 1.e-14 # tolerance sur les raccords entre courbes C1
    eps = eps*eps # car les distances sont au carres dans la suite

    # curves est maintenant une liste de i-arrays representant une segmentation
    # du profil d'un corps (continus)
    ne = len(curves)

    # Calcul des relations entre courbes C1
    ext = N.ones((ne,2), dtype=E_NpyInt)# cas d'extension
    nghb = N.zeros((ne,2), dtype=E_NpyInt); nghb[:] = -1 # no de la courbe voisine
    nghbind = N.zeros((ne,2), dtype=E_NpyInt) # index corresp. sur la courbe voisine
    # coord des pts de contact
    xpts = N.zeros((ne,2), dtype=N.float64)
    ypts = N.zeros((ne,2), dtype=N.float64)
    zpts = N.zeros((ne,2), dtype=N.float64)
    # courbure au pts de contact
    cpts = N.zeros((ne,2), dtype=N.float64)
    # normale aux points de contact
    surfxpts = N.zeros((ne,2), dtype=N.float64)
    surfypts = N.zeros((ne,2), dtype=N.float64)
    surfzpts = N.zeros((ne,2), dtype=N.float64)

    # Calcul des coord. des Pts de contact entre courbes C1
    for c in range(ne):
        s = curves[c]
        faces = P.exteriorFaces(s)
        xpts[c,0] = faces[1][0][0]
        ypts[c,0] = faces[1][1][0]
        zpts[c,0] = faces[1][2][0]
        xpts[c,1] = faces[1][0][1]
        ypts[c,1] = faces[1][1][1]
        zpts[c,1] = faces[1][2][1]
    #for c in xrange(ne):
    #    print 'courbe ',c,' x0=',xpts[c,0],' ->',xpts[c,1]

    # No des courbes voisines
    for c in range(ne):
        d = c+1
        while (d < ne):
            distx = xpts[c,0] - xpts[d,0]
            disty = ypts[c,0] - ypts[d,0]
            distz = zpts[c,0] - zpts[d,0]
            if (distx*distx+disty*disty+distz*distz < eps):
                nghb[c,0] = d
                nghbind[c,0] = 1
                nghb[d,0] = c
                nghbind[d,0] = 1
            distx = xpts[c,0] - xpts[d,1]
            disty = ypts[c,0] - ypts[d,1]
            distz = zpts[c,0] - zpts[d,1]
            if (distx*distx+disty*disty+distz*distz < eps):
                nghb[c,0] = d
                nghbind[c,0] = curves[d][2]
                nghb[d,1] = c
                nghbind[d,1] = 1
            distx = xpts[c,1] - xpts[d,0]
            disty = ypts[c,1] - ypts[d,0]
            distz = zpts[c,1] - zpts[d,0]
            if (distx*distx+disty*disty+distz*distz < eps):
                nghb[c,1] = d
                nghbind[c,1] = 1
                nghb[d,0] = c
                nghbind[d,0] = curves[c][2]
            distx = xpts[c,1] - xpts[d,1]
            disty = ypts[c,1] - ypts[d,1]
            distz = zpts[c,1] - zpts[d,1]
            if (distx*distx+disty*disty+distz*distz < eps):
                nghb[c,1] = d
                nghbind[c,1] = curves[d][2]
                nghb[d,1] = c
                nghbind[d,1] = curves[c][2]
            d += 1

    #for c in xrange(ne):
    #    print 'courbe ',c,' voisin0=',nghb[c,0],', voisin1=',nghb[c,1]

    # Courbe unique
    unique = T.join(curves)

    # Courbure sur unique
    curvature = D.getCurvatureAngle(unique)
    out = C.addVars([unique, curvature])
    #C.convertArrays2File([out], 'angle.plt')
    #unique = C.convertArray2Tetra(unique)

    # Normale sur unique
    unique2 = T.addkplane(unique)
    nunique = G.getNormalMap(unique2)
    nunique = C.center2Node(nunique)
    nunique = C.normalize(nunique, ['sx','sy','sz'])

    # On determine l'angle de courbure et les normales pour les pts de contact
    coord = unique[1]
    npts = coord.shape[1]
    for i in range(ne):
        for j in range(2):
            x = xpts[i,j]; y = ypts[i,j]; z = zpts[i,j]
            for ind in range(npts):
                distx = coord[0,ind] - x
                disty = coord[1,ind] - y
                distz = coord[2,ind] - z
                if (distx*distx+disty*disty+distz*distz < eps):
                    cpts[i,j] = curvature[1][0,ind]
                    surfxpts[i,j] = nunique[1][0,ind]
                    surfypts[i,j] = nunique[1][1,ind]
                    surfzpts[i,j] = nunique[1][2,ind]
                    break

    #for c in range(ne):
    #    print 'courbe %d: angle0=%f, angle1=%f'%(c, cpts[c,0], cpts[c,1])

    #for c in xrange(ne):
    #    print 'courbe ',c,' n0=',surfxpts[c,0],surfypts[c,0],surfzpts[c,0],\
    #          ', n1=',surfxpts[c,1],surfypts[c,1],surfzpts[c,1]

    # Calcul de l'extension
    # ext=1 ; extension chimere
    # ext=0 ; TFI paroi
    # ext=-1; TFI MD + TTM
    # ext=2; extension coincidente
    for c in range(ne):
        if (cpts[c,0] >= alphaMax): # extension coincidente
            ext[c,0] = 2
        elif (cpts[c,0] >= 185): # extension chimere
            ext[c,0] = 1
        elif (cpts[c,0] < 90+38): # TFI 2 parois
            ext[c,0] = 0
        elif (cpts[c,0] >= 90+38): # TFI MD + TTM
            ext[c,0] = -1

        if (cpts[c,1] >= alphaMax): # extension coincidente
            ext[c,1] = 2
        elif (cpts[c,1] >= 185): # extension chimere
            ext[c,1] = 1
        elif (cpts[c,1] < 90+38): # TFI 2 parois
            ext[c,1] = 0
        elif (cpts[c,1] >= 90+38): # TFI MD + TTM
            ext[c,1] = -1

    for c in range(ne):
        print ('courbe %d: ext0=%d, ext1=%d'%(c, ext[c,0], ext[c,1]))

    # Detection de la hauteur maximum admissible en fonction de la longueur
    # de chaque segment
    n = 0
    for i in curves:
        l = D.getLength(i)
        if (h > 0.9*l and ext[n,0] != -1 and ext[n,1] != -1):
            h = 0.9*l
            print("Warning: height changed to", h,"...")
            print("...because length of segment", n, "is", l)
        n = n+1

    # Detection de la hauteur maximum admissible : h < Rc, rayon de courbure
    n = 0
    for i in curves:
        r = D.getCurvatureRadius(i)
        rp = r[1]
        rmin = 1.e8
        for x in range(r[2]):
            if (rp[0,x] <= 0 and -rp[0,x] < rmin):
                rmin = -rp[0,x]

        if ( h > rmin ):
            h = float(0.99 * rmin)
            print("Warning: height changed to", h, "...")
            print("...because curvature radius in segment", n, "is", rmin)
        n = n+1

    # Detection de la densite minimum
    nj = int(h*density)+1
    if (nj < 4):
        density = 4./h
        print("Warning: density changed to", density)
    for i in curves:
        l = D.getLength(i)
        ni = int(l*density)+1
        if (ni < 2):
            density = 2./l
            print("Warning: density changed to", density)

    # Calcul automatique de l'extension
    extension = max(int(h*density)+1, 5)

    # Distribution suivant la hauteur
    nj = int(h*density)+1
    # Distribution reguliere utilisee dans la generation
    distrib = G.cart((0,0,0), (1./nj,1,1), (nj+1,1,1))
    add = max(nj * h / (20 * yplus), 1)
    add = min(int(add), 3*nj)
    add = int(addFactor*add)
    # Distribution avec resserrement utilisee a la fin pour remailler en j
    distrib2 = G.enforcePlusX(distrib, yplus/h, nj-1, add)
    nj = distrib[2]; delta = C.array('d', nj, 1, 1)
    for i in range(nj): delta[1][0,i] = h*distrib[1][0,i]
    # Generation des maillages
    mesh = []; walls = []
    cm = 0

    for c in range(ne):
        if (ext[c,0] == 1 and ext[c,1] == 1):
            m = generateExtExt(curves[c], density, extension, delta)
            m = G.map(m, distrib2, 2)
            m = T.reorder(m, (-1,2,3))
            m = T.addkplane(m)
            mesh.append(m)

            buildBC(m, walls, ext[c,0], ext[c,1], extension)
        elif (ext[c,0] == 2 and ext[c,1] == 2):
            if ne == 1: # une seule courbe
                nghb[c,0] = 0; nghb[c,1] = 0
                nghbind[c,0] = curves[c][2]; nghbind[c,1] = 1 
            m = generateExtExtMatch(c, curves, density, extension, delta,
                                    nghb[c,0], nghb[c,1],
                                    nghbind[c,0], nghbind[c,1])
            m = G.map(m, distrib2, 2)
            m = T.reorder(m, (-1,2,3))
            m = T.addkplane(m)
            mesh.append(m)
            buildBC(m, walls, ext[c,0], ext[c,1], extension)

        else:
            m = generateOtherCases(curves[c], curves,
                                   distrib, density, extension, delta, h,
                                   ext[c,0], ext[c,1],
                                   nghb[c,0], nghb[c,1],
                                   nghbind[c,0], nghbind[c,1],
                                   surfxpts[c,0], surfypts[c,0], surfzpts[c,0],
                                   surfxpts[c,1], surfypts[c,1], surfzpts[c,1],
                                   cpts[c,0], cpts[c,1])
            m = G.map(m, distrib2, 2)
            m = T.addkplane(m)
            mesh.append(m)
            buildBC(m, walls, ext[c,0], ext[c,1], extension)

    # Si une seule grille a ete engendre, on la coupe en 2
    if len(mesh) == 1:
        m = mesh[0]; isp = m[2]//2
        m1 = T.subzone(m, (1,1,1), (isp,m[3],m[4]))
        m2 = T.subzone(m, (isp,1,1), (m[2],m[3],m[4]))
        del mesh[0]; mesh.append(m1); mesh.append(m2)
        w = walls[0]; l = len(w); wf1 = []; wf2 = []
        if l > 0:
            r = w[0]
            wf1.append([r[0],isp,r[2],r[3],r[4],r[5]])
            wf2.append([1,r[1]-isp+1,r[2],r[3],r[4],r[5]])
        if l > 1:
            r = w[1]
            if r[0] == 1:
                wf1.append([r[0],r[1],r[2],r[3],r[4],r[5]])
            else:
                wf2.append([r[0]-isp+1,r[1]-isp+1,r[2],r[3],r[4],r[5]])
        if l > 2:
            r = w[2]
            if r[0] == 1:
                wf1.append([r[0],r[1],r[2],r[3],r[4],r[5]])
            else:
                wf2.append([r[0]-isp+1,r[1]-isp+1,r[2],r[3],r[4],r[5]])
        del walls[0]; walls.append(wf1); walls.append(wf2)
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
            if (ind1 == index): return ind2
            if (ind2 == index): return ind1
    return -1

#==============================================================================
# Extension lineaire
def generateExtExt(curve, density, extension, delta):

    l = D.getLength(curve)
    ni = int(l*density)+1
    hm = l/(ni-1)
    ext = extension
    vars = curve[0]
    curve[0] = 'x,y,z'
    # extension lineaire
    x1 = curve[1][0,0]; y1 = curve[1][1,0]; z1 = curve[1][2,0]
    x2 = curve[1][0,1]; y2 = curve[1][1,1]; z2 = curve[1][2,1]
    norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
    px1 = x1 - ext*(x2-x1)*hm/norm
    py1 = y1 - ext*(y2-y1)*hm/norm
    pz1 = z1 - ext*(z2-z1)*hm/norm
    line = D.line((x1,y1,z1), (px1,py1,pz1),ext)
    curve2 = T.join(line, curve)
    n = curve[1].shape[1]
    x1 = curve[1][0,n-2]; y1 = curve[1][1,n-2]; z1 = curve[1][2,n-2]
    x2 = curve[1][0,n-1]; y2 = curve[1][1,n-1]; z2 = curve[1][2,n-1]
    norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
    px1 = x2 + ext*(x2-x1)*hm/norm
    py1 = y2 + ext*(y2-y1)*hm/norm
    pz1 = z2 + ext*(z2-z1)*hm/norm
    line = D.line((x2,y2,z2),(px1,py1,pz1),ext)
    curve2 = T.join(curve2, line)

    l = D.getLength(curve2)
    ni = ni + 2*ext
    hd = 1./(ni-1)
    db = G.cart( (0,0,0), (hd,1,1), (ni,1,1) )
    m = G.map(curve2, db)
    m = G.addNormalLayers(m, delta,niter=50,check=0)
    m = T.reorder(m, (1,-3,2))
    curve[0] = vars
    return m

#==============================================================================
# Extension coincidente des deux cotes
#==============================================================================
def generateExtExtMatch(c, curves, density, extension, delta,
                        nghb1, nghb2, nghbind1, nghbind2):
    curve = curves[c]
    f = curve[1]; np = f.shape[1]
    l = D.getLength(curve)
    ni = int(l*density)+1
    hm = l/(ni-1)
    ext = extension
    vars = curve[0]
    curve[0] = 'x,y,z'
    r0 = G.cart((0,0,0), (1./(ni-1),1,1), (ni,1,1))
    curve2 = G.map(curve, r0)
    f2 = curve2[1]
    curvev1 = curves[int(nghb1)]
    lv1 = D.getLength(curvev1)
    niv1 = int(lv1*density)+1
    hmv1 = lv1/(niv1-1)

    curvev2 = curves[int(nghb2)]
    lv2 = D.getLength(curvev2)
    niv2 = int(lv2*density)+1
    hmv2 = lv2/(niv2-1)

    if nghbind1 == 1: indv1 = 1
    else: indv1 = nghbind1-2
    if nghbind2 == 1: indv2 = 1
    else: indv2 = nghbind2-2

    # extension coincidente en i=0
    x1 = f[0,0]; y1 = f[1,0]; z1 = f[2,0]
    x2 = f[0,1]; y2 = f[1,1]; z2 = f[2,1]
    norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
    px1 = x1 - ext*(x2-x1)*hm/norm
    py1 = y1 - ext*(y2-y1)*hm/norm
    pz1 = z1 - ext*(z2-z1)*hm/norm

    x3 = curvev1[1][0,indv1]; y3 = curvev1[1][1,indv1]; z3 = curvev1[1][2,indv1]
    # recup des coordonnees du pt suivant de la courbe voisine
    norm = math.sqrt( (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3) )
    px2 = x1 - ext*(x3-x1)*hmv1/norm
    py2 = y1 - ext*(y3-y1)*hmv1/norm
    pz2 = z1 - ext*(z3-z1)*hmv1/norm
    line1 = D.line((x1,y1,z1), (0.5*(px1+px2),0.5*(py1+py2),0.5*(pz1+pz2)),ext)
    curve2 = T.join(line1, curve2)

    # extension coincidente en i=n-1
    x1 = f[0,np-2]; y1 = f[1,np-2]; z1 = f[2,np-2]
    x2 = f[0,np-1]; y2 = f[1,np-1]; z2 = f[2,np-1]
    norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
    px1 = x2 + ext*(x2-x1)*hm/norm
    py1 = y2 + ext*(y2-y1)*hm/norm
    pz1 = z2 + ext*(z2-z1)*hm/norm
    x3 = curvev2[1][0,indv2]; y3 = curvev2[1][1,indv2]; z3 = curvev2[1][2,indv2]# recup des coordonnees du pt suivant de la courbe voisine
    norm = math.sqrt( (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2) )
    px2 = x2 - ext*(x3-x2)*hmv2/norm
    py2 = y2 - ext*(y3-y2)*hmv2/norm
    pz2 = z2 - ext*(z3-z2)*hmv2/norm
    line2 = D.line((x2,y2,z2), (0.5*(px1+px2),0.5*(py1+py2),0.5*(pz1+pz2)),ext)

    if nghb1 != c and nghb2 != c:
        curve2 = T.join(curve2, line2)
        m = G.addNormalLayers(curve2, delta, niter=50, check=0)
        m = T.reorder(m, (1,-3,2))
        return m
    else: # cas ou la courbe est coincidente sur elle meme
        curve = G.map(curve, r0)
        l = D.getLength(curve)
        curve2 = T.subzone(curve, (2,1,1), (curve[2]-1,1,1))
        line1 = T.join(line1, T.subzone(curve,(1,1,1),(2,1,1)))
        line2 = T.join(T.subzone(curve,(curve[2]-1,1,1),(curve[2],1,1)),line2)
        m = G.addNormalLayers(curve2, delta, niter=50, check=0)
        m = T.reorder(m, (1,-3,2))
        # construction de l'extension en i=1
        d1 = line1
        d4 = T.subzone(m,(1,1,1),(1,m[3],1))
        m1 = G.addNormalLayers(d1, delta, niter=50, check=0)
        m1 = T.reorder(m1, (1,-3,2))
        d3 = T.subzone(m1,(1,1,1),(1,m1[2],1))
        d3 = T.reorder(d3,(-2,1,3))
        d4 = T.reorder(d4,(-2,1,3))
        ind = d3[2]-1
        px1 = d3[1][0,ind]; py1 = d3[1][1,ind]; pz1 = d3[1][2,ind]
        px2 = d4[1][0,ind]; py2 = d4[1][1,ind]; pz2 = d4[1][2,ind]
        d2 = D.line((px1,py1,pz1),(px2,py2,pz2),d1[2])
        m1 = G.TFI([d1,d2,d3,d4])
        # construction de l'extension en i=imax
        d1 = line2
        d4 = T.subzone(m,(m[2],1,1),(m[2],m[3],1))
        m2 = G.addNormalLayers(d1, delta, niter=50, check=0)
        m2 = T.reorder(m2, (1,-3,2))
        d3 = T.subzone(m2,(m2[2],1,1),(m2[2],m2[2],1))
        d3 = T.reorder(d3,(-2,1,3))
        d4 = T.reorder(d4,(-2,1,3))
        ind = d3[2]-1
        px1 = d3[1][0,ind]; py1 = d3[1][1,ind]; pz1 = d3[1][2,ind]
        px2 = d4[1][0,ind]; py2 = d4[1][1,ind]; pz2 = d4[1][2,ind]
        d2 = D.line((px1,py1,pz1),(px2,py2,pz2),d1[2])
        m2 = G.TFI([d1,d2,d3,d4])
        m = T.join(m1,m); m = T.join(m,m2)
        return m

#==============================================================================
# Autres cas d'extension
def generateOtherCases(curve, curves,
                       distrib, density, extension, delta, h,
                       ext1, ext2, nghb1, nghb2, nghbind1, nghbind2,
                       n0x, n0y, n0z, n1x, n1y, n1z, angle1, angle2):
    vars = curve[0]
    curve2 = T.addkplane(curve)
    n = G.getNormalMap(curve2)
    n = C.center2Node(n)
    n = C.normalize(n, ['sx','sy','sz'])
    f = curve[1]
    l = D.getLength(curve)
    ni = int(l*density)+1
    hm = l/(ni-1);
    ext = extension
    r0 = G.cart((0,0,0), (1./(ni-1),1,1), (ni,1,1))
    curve0 = G.map(curve,r0)
    f0 = curve0[1]
    curve[0] = 'x,y,z'
    d2 = []
    if (ext1 == 2):
        curvev1 = curves[int(nghb1)]
        if nghbind1 == 1: indv1 = 1
        else: indv1 = nghbind1-2
        lv1 = D.getLength(curvev1)
        niv1 = int(lv1*density)+1
        hmv1 = lv1/(niv1-1)

        curvev2 = curves[int(nghb2)]
        if nghbind2 == 1: indv2 = 1
        else: indv2 = nghbind2-2

        x1 = f[0,0]; y1 = f[1,0]; z1 = f[2,0]
        x2 = f[0,1]; y2 = f[1,1]; z2 = f[2,1]
        d1 = curve0; nid1 = d1[2]
        px3 = x1 + n[1][0,0] * h
        py3 = y1 + n[1][1,0] * h
        pz3 = z1 + n[1][2,0] * h
        d3 = D.line((x1,y1,z1), (px3,py3,pz3),ext)

        if (ext2 == 0):# TFI paroi            
            if (nghbind2 == 1):
                h0 = h / math.sin( angle2 * math.pi / 180. )
                h0 = min(h0, 2*h)
                ind = D.getDistantIndex(curvev2, 1, h0)
                d4 = T.subzone( curvev2, (1,1,1), (ind,1,1) )
            else:
                h0 = h / math.sin( angle2 * math.pi / 180. )
                h0 = min(h0, 2*h)
                ind = D.getDistantIndex(curvev2, nghbind2, -h0)
                d4 = T.subzone( curvev2, (ind,1,1), (curvev[2],1,1) )
                d4 = T.reorder(d4, (-1,2,3))

        elif (ext2 == -1):# TFI MD + TTM
            nil = d1[2]
            px2 = d1[1][0,nil-1]
            py2 = d1[1][1,nil-1]
            pz2 = d1[1][2,nil-1]
            px4 = px2 + n1x * h
            py4 = py2 + n1y * h
            pz4 = pz2 + n1z * h
            d4 = D.line((px2,py2,pz2), (px4,py4,pz4),ext)

        elif (ext2 == 1): #extension recouvrante
            p = f.shape[1]
            x1 = f[0,p-1]; y1 = f[1,p-1]; z1 = f[2,p-1]
            x2 = f[0,p-2]; y2 = f[1,p-2]; z2 = f[2,p-2]
            norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
            px1 = x1 + ext*(x1-x2)*hm/norm
            py1 = y1 + ext*(y1-y2)*hm/norm
            pz1 = z1 + ext*(z1-z2)*hm/norm
            line = D.line((x1,y1,z1), (px1,py1,pz1),ext)
            d1 = T.join(d1, line)
            px3 = px1 + n[1][0,p-1] * h
            py3 = py1 + n[1][1,p-1] * h
            pz3 = pz1 + n[1][2,p-1] * h
            d4 = D.line((px1,py1,pz1), (px3,py3,pz3),ext)

        d2 = buildd2(d1, d3, d4, h, ext1, ext2, angle1, angle2) 
        r02 = G.cart((0,0,0), (1./(d1[2]-1),1,1), (d1[2],1,1))
        d2 = G.map(d2,r02)

        # Ajout de l'extension coincidente en i=1
        x1 = f[0,0]; y1 = f[1,0]; z1 = f[2,0]
        x2 = f[0,1]; y2 = f[1,1]; z2 = f[2,1]
        norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
        px1 = x1 - ext*(x2-x1)*hm/norm
        py1 = y1 - ext*(y2-y1)*hm/norm
        pz1 = z1 - ext*(z2-z1)*hm/norm
        x3 = curvev1[1][0,indv1]; y3 = curvev1[1][1,indv1]; z3 = curvev1[1][2,indv1]# recup des coordonnees du pt suivant de la courbe voisine
        norm = math.sqrt( (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3) )
        px2 = x1 - ext*(x3-x1)*hmv1/norm
        py2 = y1 - ext*(y3-y1)*hmv1/norm
        pz2 = z1 - ext*(z3-z1)*hmv1/norm
        line1 = D.line((x1,y1,z1), (0.5*(px1+px2),0.5*(py1+py2),0.5*(pz1+pz2)),ext)        
        line1 = T.reorder(line1,(-1,2,3))
        # ajout du nouveau d3 pres de l'extension coincidente
        line4 = d3
        m1 = G.addNormalLayers(line1, delta,niter=50, check=0)
        m1 = T.reorder(m1, (1,-3,2))
        d3 = T.subzone(m1,(1,1,1),(1,m1[2],1)); d3 = T.reorder(d3,(-2,1,3))
        ind = line4[2]-1
        px1 = line4[1][0,ind]; py1 = line4[1][1,ind]; pz1 = line4[1][2,ind]
        px2 = d3[1][0,ind]; py2 = d3[1][1,ind]; pz2 = d3[1][2,ind]
        line2 = D.line((px1,py1,pz1),(px2,py2,pz2),ext)
        d1 = T.join(line1,d1)
        d2 = T.join(line2,d2)

    elif (ext1 == 1):
        x1 = f[0,0]; y1 = f[1,0]; z1 = f[2,0]
        x2 = f[0,1]; y2 = f[1,1]; z2 = f[2,1]
        norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
        px1 = x1 - ext*(x2-x1)*hm/norm
        py1 = y1 - ext*(y2-y1)*hm/norm
        pz1 = z1 - ext*(z2-z1)*hm/norm
        line = D.line((x1,y1,z1), (px1,py1,pz1),ext)
        d1 = T.join(line, curve0)

        px3 = px1 + n[1][0,0] * h
        py3 = py1 + n[1][1,0] * h
        pz3 = pz1 + n[1][2,0] * h
        d3 = D.line((px1,py1,pz1), (px3,py3,pz3),ext)
        #print 'pts',px1,py1,pz1,'->',px3,py3,pz3

        if (ext2 == 0):
            curvev = curves[int(nghb2)]
            if (nghbind2 == 1):
                h0 = h / math.sin( angle2 * math.pi / 180. )
                h0 = min(h0, 2*h)
                ind = D.getDistantIndex(curvev, 1, h0)
                d4 = T.subzone( curvev, (1,1,1), (ind,1,1) )
            else:
                h0 = h / math.sin( angle2 * math.pi / 180. )
                h0 = min(h0, 2*h)
                ind = D.getDistantIndex(curvev, nghbind2, -h0)
                d4 = T.subzone( curvev, (ind,1,1), (curvev[2],1,1) )
                d4 = T.reorder(d4, (-1,2,3))

        elif (ext2 == -1):
            nil = d1[2]
            px2 = d1[1][0,nil-1]
            py2 = d1[1][1,nil-1]
            pz2 = d1[1][2,nil-1]
            px4 = px2 + n1x * h
            py4 = py2 + n1y * h
            pz4 = pz2 + n1z * h
            d4 = D.line((px2,py2,pz2), (px4,py4,pz4),ext)

        elif (ext2 == 2):
            d1,d4, d2 = buildd1d4Match(curves[int(nghb2)], nghbind2, density, delta, \
                                       f, ext, d1, d3, hm, h, n, ext1, ext2, angle1, angle2)

    elif (ext1 == 0):
        curvev = curves[int(nghb1)]
        if (nghbind1 == 1):
            h0 = comph0(h, angle1)
            ind = D.getDistantIndex(curvev, 1, h0)
            d3 = T.subzone( curvev, (1,1,1), (ind,1,1) )
        else:
            h0 = comph0(h, angle1)
            ind = D.getDistantIndex(curvev, nghbind1, -h0)
            if ( ind > 0 ):
                nil = curvev[2]
                d3 = T.subzone( curvev, (ind,1,1), (nil,1,1) )
            else:
                d3 = C.copy(curvev)
            d3 = T.reorder(d3, (-1,2,3))

        if (ext2 == 1):
            p = f.shape[1]
            x1 = f[0,p-1]; y1 = f[1,p-1]; z1 = f[2,p-1]
            x2 = f[0,p-2]; y2 = f[1,p-2]; z2 = f[2,p-2]
            norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
            px1 = x1 + ext*(x1-x2)*hm/norm
            py1 = y1 + ext*(y1-y2)*hm/norm
            pz1 = z1 + ext*(z1-z2)*hm/norm
            line = D.line((x1,y1,z1), (px1,py1,pz1),ext)
            d1 = T.join(curve0, line)
            px3 = px1 + n[1][0,p-1] * h
            py3 = py1 + n[1][1,p-1] * h
            pz3 = pz1 + n[1][2,p-1] * h
            d4 = D.line((px1,py1,pz1), (px3,py3,pz3),ext)

        elif (ext2 == 2):
            d1,d4, d2 = buildd1d4Match(curves[int(nghb2)], nghbind2, density, delta, \
                                       f, ext, d1, d3, hm, h, n, ext1, ext2, angle1, angle2)

        elif (ext2 == 0):
            d1 = curve0
            curvev = curves[int(nghb2)]
            if (nghbind2 == 1):
                h0 = comph0(h, angle2)
                ind = D.getDistantIndex(curvev, 1, h0)
                d4 = T.subzone( curvev, (1,1,1), (ind,1,1) )
            else:
                h0 = comph0(h, angle2)
                ind = D.getDistantIndex(curvev, nghbind2, -h0)
                nil = curvev[2]
                d4 = T.subzone( curvev, (ind,1,1), (nil,1,1) )
                d4 = T.reorder(d4, (-1,2,3))

        else: # ext2 = -1
            d1 = curve0
            p = f.shape[1]
            x1 = f[0,p-1]; y1 = f[1,p-1]; z1 = f[2,p-1]
            px3 = x1 + n1x * h
            py3 = y1 + n1y * h
            pz3 = z1 + n1z * h
            d4 = D.line((x1,y1,z1), (px3,py3,pz3),ext)

    else: # ext1 == -1
        d1 = curve0
        px1 = d1[1][0,0]
        py1 = d1[1][1,0]
        pz1 = d1[1][2,0]
        px3 = px1 + n0x * h
        py3 = py1 + n0y * h
        pz3 = pz1 + n0z * h
        #print 'points ',px1,py1,pz1,'->',px3,py3,pz3
        d3 = D.line((px1,py1,pz1), (px3,py3,pz3),ext)

        if (ext2  == 1):
            nil = f.shape[1]
            x1 = f[0,nil-1]; y1 = f[1,nil-1]; z1 = f[2,nil-1]
            x2 = f[0,nil-2]; y2 = f[1,nil-2]; z2 = f[2,nil-2]
            norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) )
            px1 = x1 + ext*(x1-x2)*hm/norm
            py1 = y1 + ext*(y1-y2)*hm/norm
            pz1 = z1 + ext*(z1-z2)*hm/norm
            line = D.line((x1,y1,z1), (px1,py1,pz1),ext)
            d1 = T.join(d1, line)
            px3 = px1 + n[1][0,nil-1] * h
            py3 = py1 + n[1][1,nil-1] * h
            pz3 = pz1 + n[1][2,nil-1] * h
            d4 = D.line((px1,py1,pz1), (px3,py3,pz3),ext)

        elif (ext2 == 2):
            d1,d4,d2 = buildd1d4Match(curves[int(nghb2)], nghbind2, density, delta, \
                                      f, ext, d1, d3, hm, h, n, ext1, ext2, angle1, angle2)
        if (ext2 == 0):
            curvev = curves[int(nghb2)]
            if (nghbind2 == 1):
                h0 = comph0(h, angle2)
                ind = D.getDistantIndex(curvev, 1, h0)
                #print 'indice lointain ', ind
                d4 = T.subzone( curvev, (1,1,1), (ind,1,1) )
            else:
                h0 = comph0(h, angle2)
                ind = D.getDistantIndex(curvev, nghbind2, -h0)
                nil = curvev[2]
                d4 = T.subzone( curvev, (ind,1,1), (nil,1,1) )
                d4 = T.reorder(d4, (-1,2,3))

        if (ext2 == -1):
            nil = d1[2]
            px2 = d1[1][0,nil-1]
            py2 = d1[1][1,nil-1]
            pz2 = d1[1][2,nil-1]
            px4 = px2 + n1x * h
            py4 = py2 + n1y * h
            pz4 = pz2 + n1z * h
            d4 = D.line((px2,py2,pz2), (px4,py4,pz4),ext)

    vars0 = 'x,y,z'
    d1[0] = vars0; d3[0] = vars0; d4[0] = vars0
    if d2 == []: d2 = buildd2(d1, d3, d4, h, ext1, ext2, angle1, angle2)
    d1[0] = vars; d2[0] = vars; d3[0] = vars; d4[0] = vars; distrib[0] = vars

    if ext1 > 0: ni = ni+ext
    if ext2 > 0: ni = ni+ext

    r1 = d1; r2 = d2
    r3 = G.map(d3, distrib)
    r4 = G.map(d4, distrib)
    #r1 = T.reorder(r1, (-1,2,3))
    #r3 = T.reorder(r3, (-1,2,-3))
    m = G.TFI([r3, r4, r1, r2])
    m = T.reorder(m, (-1,2,3))

    #m = G.TTM(m, 20)   
    curve[0] = vars
    return m

#=============================================================================
# Construction de d4 et modification de d1 dans le cas d extension coincidente
#=============================================================================
def buildd1d4Match(curvev2, nghbind2, density, delta, f, ext, d1, d3, hm, h, n, ext1, ext2, angle1, angle2):
    if nghbind2 == 1: indv2 = 1
    else: indv2 = nghbind2-2
    lv2 = D.getLength(curvev2)
    niv2 = int(lv2*density)+1
    hmv2 = lv2/(niv2-1)
    nid1 = d1[2]
    r0 = G.cart((0,0,0), (1./(nid1-1),1,1), (nid1,1,1))
    p = f.shape[1]
    x1 = f[0,p-2]; y1 = f[1,p-2]; z1 = f[2,p-2]
    x2 = f[0,p-1]; y2 = f[1,p-1]; z2 = f[2,p-1]

    px = x2 + n[1][0,p-1] * h
    py = y2 + n[1][1,p-1] * h
    pz = z2 + n[1][2,p-1] * h
    d4 = D.line((x2,y2,z2), (px,py,pz),ext)
    d2 = buildd2(d1, d3, d4, h, ext1, ext2, angle1, angle2)
    d2 = G.map(d2,r0)
    #
    # recup des coordonnees du pt suivant de la courbe voisine
    #
    norm = math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
    px1 = x2 + ext*(x2-x1)*hm/norm
    py1 = y2 + ext*(y2-y1)*hm/norm
    pz1 = z2 + ext*(z2-z1)*hm/norm
    #
    x3 = curvev2[1][0,indv2]; y3 = curvev2[1][1,indv2]; z3 = curvev2[1][2,indv2]
    norm3 = math.sqrt( (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2) )
    #
    px2 = x2 + ext*(x2-x3)*hmv2/norm3
    py2 = y2 + ext*(y2-y3)*hmv2/norm3
    pz2 = z2 + ext*(z2-z3)*hmv2/norm3

    px1 = 0.5*(px1+px2); py1 = 0.5*(py1+py2); z1 = 0.5*(pz1+pz2)
    line1 = D.line((x2,y2,z2), (px1,py1,pz1),ext)
    #
    m1 = G.addNormalLayers(line1, delta,niter=50, check=0)
    m1 = T.reorder(m1, (1,-3,2))
    line4 = T.subzone(m1,(m1[2],1,1),(m1[2],m1[2],1),ext)
    line4 = T.reorder(line4,(-2,1,3))
    line3 = d4
    ind = line3[2]-1
    px1 = line3[1][0,ind]; py1 = line3[1][1,ind]; pz1 = line3[1][2,ind]
    px2 = line4[1][0,ind]; py2 = line4[1][1,ind]; pz2 = line4[1][2,ind]
    line2 = D.line((px1,py1,pz1),(px2,py2,pz2),ext)
    d1 = T.join(d1,line1)
    d2 = T.join(d2,line2)
    d4 = line4
    return d1,d4,d2

#==============================================================================
def buildd2(d1, d3, d4, h, ext1, ext2, angle1, angle2):

    # Abscisse curviligne sur d1
    s = D.getCurvilinearAbscissa(d1)

    # normale a d1
    d1k = T.addkplane(d1)
    n = G.getSmoothNormalMap(d1k, niter=100)
    n = C.center2Node(n)
    n = C.normalize(n, ['sx','sy','sz'])
    np = n[1]

    # Amortissement
    l = D.getLength(d1)
    if (ext1 == 0):
        H = compH(h, angle1)
        s1 = H / l
        s1 = min(s1, 0.5)
    elif (ext1 >0):
        s1 = 1./float(d1[2])
    elif (ext1 == -1 and angle1 >= 180.):
        s1 = 1./float(d1[2])
    else:
        s1 = 1./18.

    if (ext2 == 0):
        H = compH(h, angle2)
        s2 = 1.- H / l
        s2 = max(s2, 0.5)
    elif (ext2 >0):
        s2 = (d1[2]-1.)/float(d1[2])
    elif (ext2 == -1 and angle2 >= 180.):
        s2 = (d1[2]-1.)/float(d1[2])
    else:
        s2 = 17./18.

    ind13 = D.getDistantIndex(d1, 1, s1*l)
    ind23 = D.getDistantIndex(d1, 1, s2*l)
    #print ext1,ext2,s1,s2, ind13, ind23, d1[2]
    #print d1[1][0,0], d1[1][1,0] 
    #print ind13, ind23, l, 1./3.*l

    d1p = d1[1]; d3p = d3[1]; d4p = d4[1]

    x1 = d1p[0,0]; y1 = d1p[1,0]; z1 = d1p[2,0]
    p = d3[2]
    x3 = d3p[0,p-1]; y3 = d3p[1,p-1]; z3 = d3p[2,p-1]
    p = d1[2]
    x2 = d1p[0,p-1]; y2 = d1p[1,p-1]; z2 = d1p[2,p-1]
    p = d4[2]
    x4 = d4p[0,p-1]; y4 = d4p[1,p-1]; z4 = d4p[2,p-1]

    h1x = x3 - x1; h1y = y3 - y1; h1z = z3 - z1
    h2x = x4 - x2; h2y = y4 - y2; h2z = z4 - z2
    hmx = 0.5*(h1x+h2x); hmy = 0.5*(h1y+h2y); hmz = 0.5*(h1z+h2z)

    ni = s[2]; sp = s[1]
    vect = C.array('hx,hy,hz', ni, 1, 1)

    for i in range(ind13):
        sl = 1./s1*sp[0,i]
        vect[1][0,i] = (1.-sl)*h1x + sl*h*np[0,ind13]
        vect[1][1,i] = (1.-sl)*h1y + sl*h*np[1,ind13]
        vect[1][2,i] = (1.-sl)*h1z + sl*h*np[2,ind13]
    for i in range(ind13, ind23):
        vect[1][0,i] = h*np[0,i]
        vect[1][1,i] = h*np[1,i]
        vect[1][2,i] = h*np[2,i]
    for i in range(ind23, ni):
        sl = (sp[0,i]-s2)/(1. -s2)
        vect[1][0,i] = sl*h2x + (1.-sl)*h*np[0,ind23]
        vect[1][1,i] = sl*h2y + (1.-sl)*h*np[1,ind23]
        vect[1][2,i] = sl*h2z + (1.-sl)*h*np[2,ind23]  
    d2 = T.deform(d1, vect)
    return d2

#=============================================================================
# Build boundary conditions
#=============================================================================
def buildBC(m, walls,  ext1, ext2, extension):

    # Walls
    wl = []
    if ext1 > 0 and ext2 > 0:
        if ext1 == 1 and ext2 == 1:
            i1 = extension+1; i2 = m[2]-extension
        else:
            i1 = extension; i2 = m[2]-extension+1
    elif ext1 <= 0 and ext2 <= 0:
        i1 = 1; i2 = m[2]
    elif ext1 <= 0 and ext2 > 0:
        i2 = m[2];  i1 = extension
    else:
        i2 = m[2]-extension+1; i1 = 1  

    wrange = [i1,i2,1,1,1,m[4]]
    wl.append(wrange)

    if (ext2 == 0):
        wrange = [1,1,1,m[3],1,m[4]]
        wl.append(wrange)

    if (ext1 == 0):
        wrange = [m[2],m[2],1,m[3],1,m[4]]
        wl.append(wrange)

    walls.append(wl)
    return

#==============================================================================
def comph0(h, angle):
    if (abs(angle) < 1.e-6):
        return 2*h
    if (angle <= 90 and angle >= 0):
        h0 = h / math.sin( angle * math.pi / 180. )
    else:
        h0 = h / math.cos( (angle-90.) * math.pi / 180. )
    h0 = min(h0, 2*h)
    return h0

#==============================================================================
def compH(h, angle):
    if (angle <= 90 and angle >= 0):
        H = h / math.tan( angle * math.pi / 180. )
    else:
        H = h * math.tan( (angle-90.) * math.pi / 180. )
    return 3*abs(H)
