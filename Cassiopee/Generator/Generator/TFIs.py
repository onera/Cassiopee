# Various specific TFIs
from . import Generator as G
from . import generator
try: import Converter as C
except ImportError:
    raise ImportError("TFIs: requires Converter module.")
import KCore.Vector as Vector

try: range = xrange
except: pass

import numpy

#==============================================================================
# Evalue la qualite du maillage m
# Plus le score est elevee pour le maillage est mauvais
#==============================================================================
def quality(meshes):
    """Measure the quality of mesh."""
    score = 0.
    for m in meshes:
        ortho = G.getOrthogonalityMap(m)
        vol = G.getVolumeMap(m)
        min1 = C.getMinValue(ortho, 'orthogonality')
        max1 = C.getMaxValue(ortho, 'orthogonality')
        min2 = C.getMinValue(vol, 'vol')
        max2 = C.getMaxValue(vol, 'vol')
        score = max(score, min1); score = max(score, max1)
        if min2 < 1.e-12 and max2 > 0: score += 1000.
        elif max2 < 1.e-12 and min2 > 0: score += 1000.
    return score

#==============================================================================
# Order a set of structured edges in a loop
#==============================================================================
def orderEdges(edges, tol=1.e-10):
    """Order edges in a loop."""
    import Transform
    out = [] # ordered list of edges
    pool = edges[:] # list copy
    cur = pool[0]; pool.pop(0)
    out.append(cur)
    P1p = (cur[1][0,-1],cur[1][1,-1],cur[1][2,-1])
    while len(pool) > 0:
        found = False
        for c, p in enumerate(pool):
            P0 = (p[1][0,0],p[1][1,0],p[1][2,0])
            P1 = (p[1][0,-1],p[1][1,-1],p[1][2,-1])
            if Vector.squareDist(P1p,P0) < tol*tol:
                cur = p; out.append(cur); P1p = P1; pool.pop(c); found=True; break
            if Vector.squareDist(P1p,P1) < tol*tol:
                cur = Transform.reorder(p,(-1,2,3)); out.append(cur); P1p = P0; pool.pop(c); found=True; break
        if not found: break
    return out

#==============================================================================
# Compute junction angles for already ordered edges
#==============================================================================
def computeJunctionAngles(edges):
    """Compute junction angle in radians."""
    angles = []
    for c, p in enumerate(edges):
        P0 = (p[1][0,0],p[1][1,0],p[1][2,0])
        P1 = (p[1][0,1],p[1][1,1],p[1][2,1])
        if c == 0: pn = edges[-1]
        else: pn = edges[c-1]
        P2 = (pn[1][0,-2],pn[1][1,-2],pn[1][2,-2])
        v1 = Vector.sub(P1, P0)
        v2 = Vector.sub(P2, P0)
        v1 = Vector.normalize(v1)
        v2 = Vector.normalize(v2)
        s = Vector.dot(v1, v2)
        v = Vector.norm(Vector.cross(v1, v2))
        angle = float(numpy.arctan2(v, s))
        angles.append(angle)
    return angles

#==============================================================================
# Try to merge edges to get less edges, thresold in radians, ordered edges
#==============================================================================
def mergeEdgesAngle(edges, thresold=0.3):
    """Merge edges with a junction angle less than thresold."""
    import Transform
    angles = computeJunctionAngles(edges)
    print('angles=', angles)
    medges = []
    for c, e in enumerate(edges):
        tangle = abs(angles[c] - numpy.pi)
        if tangle < thresold and len(medges) > 0:
            e2 = Transform.join(medges[-1], e); medges[-1] = e2            
        else: medges.append(e)
    tangle0 = abs(angles[0] - numpy.pi)
    if tangle0 < thresold and len(medges) > 0:
        e2 = Transform.join(medges[-1], medges[0]); medges[-1] = e2; medges.pop(0)
    print('out', len(medges))
    return medges

def findSplitInEdges(edges, N2):
    isp = 0; N = 0
    for e in edges:
        N += C.getNPts(e)-1 
        if N == N2-1: break
        if N > N2: isp = -1; break
        isp += 1
    return isp

#==============================================================================
# merge edges to get a set of 4 edges, ordered edges
# if fail, return None
#==============================================================================
def mergeEdges4(iedges):
    import Transform
    edges = iedges[:] # copie
    # total number of points
    Ntot = 0
    for e in edges: Ntot += C.getNPts(e)-1
    Ntot += 1
    N2 = Ntot//2+1
    print("Ntot=", Ntot,"N2=",N2)
    print(Ntot - (N2-1)*2+1)
    if Ntot - ((N2-1)*2+1) != 0: return None # can not be merged in 4
    
    # Find the nearest 90 degrees junction
    #angles = computeJunctionAngles(edges)
    #c90 = 0
    #for c, a in angles:
    #    if abs(a-numpy/2.) < abs(angles[c90]-numpy.pi/2): c90 = c
    #if c90 > 0: 
    #    e1 = edges[0:c90]; e2 = edges[c90:]; edges = e2+e1

    # try to find first split
    nedges = len(edges)
    roll = 0; isp = -1
    while (isp == -1 and roll < nedges):
        isp = findSplitInEdges(edges, N2)
        if isp == -1:
            edges = edges[1:]+edges[0]
            roll += 1

    print("isp=",isp,"roll=",roll,flush=True)

    # try to find second split
    e1 = edges[0:isp+1]; e2 = edges[isp+1:]
    n1 = len(e1); n2 = len(e2)
    print("resulting",n1,n2,flush=True)
    if n1 == 2: # obvious
        N1 = C.getNPts(e1[0])
        isp = findSplitInEdges(e2, N1)
        if isp == -1: print("This is a problem")
        e3 = e2[0:isp+1]; e4 = e2[isp+1:]
        e3 = Transform.join(e3); e4 = Transform.join(e4)
        return e1+[e3,e4]
    if n2 == 2: # obvious
        N2 = C.getNPts(e2[0])
        isp = findSplitInEdges(e1, N2)
        if isp == -1: print("This is a problem")
        e3 = e1[0:isp+1]; e4 = e1[isp+1:]
        e3 = Transform.join(e3); e4 = Transform.join(e4)
        return [e3,e4]+e2
    # all need to join but how?
    return None

#===============================================================================
# Build a TFI for a set of edges
# IN: edges: list of arrays defining a loop
# OUT: list of surface meshes
#===============================================================================
def allTFI__(edges):
    nedges = len(edges)
    if nedges == 4:
        try: return [TFI(edges)]
        except: pass
    elif nedges == 1:
        try: return TFIO(edges[0])
        except: pass
    elif nedges == 2:
        try: return TFIHalfO(edges[0], edges[1])
        except: pass
    elif nedges == 3:
        try: return TFITri(edges[0], edges[1], edges[2])
        except: pass
    return None

def allTFI(edges):
    print("direct TFI.")
    ret = allTFI__(edges)
    if ret is not None: return ret

    # try to merge in 4
    medges = mergeEdges4(edges)
    if medges is not None:
        print("merged4.")
        ret = allTFI__(medges)
        if ret is not None: return ret

    # Try to merge with angles
    medges = mergeEdgesAngle(edges, 0.5)
    print("merged angles ", len(medges))
    ret = allTFI__(medges)
    if ret is not None: return ret

    # Fallback use TFIStar
    #return Generator.TFIStar2(edges)
    return TFIStar(medges)

#==============================================================================
def TFI(arrays):
    """Generate a transfinite interpolation mesh from boundaries.
    Usage: TFI(arrays)"""
    return generator.TFI(arrays)

#==============================================================================
# IN: a1,a2,a3: les 3 cotes du triangle (N3-N2+N1 impair)
# OUT: 3 maillages
#==============================================================================
def TFITri(a1, a2, a3, tol=1.e-6):
    """Generate a transfinite interpolation mesh from 3 input curves."""
    import Geom as D
    import Transform as T
    N1 = a1[2]; N2 = a2[2]; N3 = a3[2]
    
    # Verif de N
    Nt = N3-N2+N1+1
    if Nt//2-Nt*0.5 != 0: raise ValueError("TFITri: N3-N2+N1 must be odd.")
    N = Nt//2
    if N < 2: raise ValueError("TFITri: invalid number of points for this operation.")
    if N > N1-1: raise ValueError("TFITri: invalid number of points for this operation.")
    if N > N2-1: raise ValueError("TFITri: invalid number of points for this operation.")

    # Assure la continuite
    P0 = (a1[1][0,N1-1], a1[1][1,N1-1], a1[1][2,N1-1])
    P00 = (a2[1][0,0], a2[1][1,0], a2[1][2,0])
    P01 = (a2[1][0,N2-1], a2[1][1,N2-1], a2[1][2,N2-1])
    if abs(P0[0]-P00[0]) + abs(P0[1]-P00[1]) + abs(P0[2]-P00[2]) > tol and abs(P0[0]-P01[0]) + abs(P0[1]-P01[1]) + abs(P0[2]-P01[2]) > tol:
        t = a2; a2 = a3; a3 = t
        N2 = a2[2]; N3 = a3[2]
        P00 = (a2[1][0,0], a2[1][1,0], a2[1][2,0])
        
    if abs(P0[0]-P00[0]) > tol: a2 = T.reorder(a2, (-1,2,3))
    elif abs(P0[1]-P00[1]) > tol: a2 = T.reorder(a2, (-1,2,3))
    elif abs(P0[2]-P00[2]) > tol: a2 = T.reorder(a2, (-1,2,3))
    P0 = (a2[1][0,N2-1], a2[1][1,N2-1], a2[1][2,N2-1])
    P00 = (a3[1][0,0], a3[1][1,0], a3[1][2,0])
    if abs(P0[0]-P00[0]) > tol: a3 = T.reorder(a3, (-1,2,3))
    elif abs(P0[1]-P00[1]) > tol: a3 = T.reorder(a3, (-1,2,3))
    elif abs(P0[2]-P00[2]) > tol: a3 = T.reorder(a3, (-1,2,3))
    #C.convertArrays2File([a1,a2,a3], 'order.plt')

    # Center
    w1 = C.array('weight', N1, 1, 1)
    w1 = C.initVars(w1, 'weight', 1)
    w2 = C.array('weight', N2, 1, 1)
    w2 = C.initVars(w2, 'weight', 1)
    w3 = C.array('weight', N3, 1, 1)
    w3 = C.initVars(w3, 'weight', 1)
    CC = G.barycenter([a1,a2,a3], [w1,w2,w3])
    
    # Subzones
    s1 = T.subzone(a1, (1,1,1), (N1-N+1,1,1))
    s2 = T.subzone(a1, (N1-N+1,1,1), (N1,1,1))
    s3 = T.subzone(a2, (1,1,1), (N2-N1+N,1,1))
    s4 = T.subzone(a2, (N2-N1+N,1,1), (N2,1,1))
    s5 = T.subzone(a3, (1,1,1), (N,1,1))
    s6 = T.subzone(a3, (N,1,1), (N3,1,1))

    # Lines
    index = N1-N
    P01 = (a1[1][0,index], a1[1][1,index], a1[1][2,index])
    index = N2-N1+N-1
    P12 = (a2[1][0,index], a2[1][1,index], a2[1][2,index])
    index = N-1
    P23 = (a3[1][0,index], a3[1][1,index], a3[1][2,index])
    l1 = D.line(CC, P01, N=N2-N1+N)
    l2 = D.line(CC, P12, N=N)
    l3 = D.line(CC, P23, N=N1-N+1)

    # TFIs
    m1 = G.TFI([s1,l1,l3,s6])
    m2 = G.TFI([l1,s2,s3,l2])
    m3 = G.TFI([l2,s4,s5,l3])

    #return [l1,l2,l3,s1,s2,s3,s4,s5,s6]
    #return [s1,l1,l3,s6,m1]
    return [m1,m2,m3]

#==============================================================================
# Cree un maillage en O avec un carre au milieu (5 maillages)
# IN: une seule courbe avec nombre impair de points
#==============================================================================
def TFIO__(a, weight, offset=0):
    import Transform as T
    import Geom as D
    
    Nt = a[2]
    
    # Calcul des points P1, P2, P3, P4
    w = C.array('weight', Nt, 1, 1)
    w = C.initVars(w, 'weight', 1.); w[1][0,0:Nt//4+1] = weight
    P1 = G.barycenter(a, w)
    w = C.initVars(w, 'weight', 1.); w[1][0,Nt//4:Nt//2+1] = weight
    P2 = G.barycenter(a, w)
    w = C.initVars(w, 'weight', 1.); w[1][0,Nt//2:3*Nt//4+1] = weight
    P3 = G.barycenter(a, w)
    w = C.initVars(w, 'weight', 1.); w[1][0,3*Nt//4:Nt] = weight
    P4 = G.barycenter(a, w)
    
    # Calcul de P'1: projete de P1 sur le cercle
    b = C.convertArray2Hexa(a); b = G.close(b)
    PP = C.array('x,y,z', 4,1,1)
    C.setValue(PP, 0, (P1[0], P1[1], P1[2]))
    C.setValue(PP, 1, (P2[0], P2[1], P2[2]))
    C.setValue(PP, 2, (P3[0], P3[1], P3[2]))
    C.setValue(PP, 3, (P4[0], P4[1], P4[2]))

    PPP = T.projectOrtho(PP, [b])
    PPP = PPP[1]
    
    PP1 = (PPP[0,0], PPP[1,0], PPP[2,0])
    PP2 = (PPP[0,1], PPP[1,1], PPP[2,1])
    
    indexPP1 = D.getNearestPointIndex(a, PP1)[0]+offset
    if indexPP1 < 0: indexPP1 = Nt+indexPP1-1
    if indexPP1 > Nt-1: indexPP1 = indexPP1-Nt-1
    
    # Renumerote a a partir de PP1
    b = T.subzone(a, (1,1,1), (indexPP1+1,1,1) )
    c = T.subzone(a, (indexPP1+1,1,1), (Nt,1,1) )
    a = T.join(c, b)
    indexPP1 = 0
    PP1 = (a[1][0,indexPP1], a[1][1,indexPP1], a[1][2,indexPP1])
    
    indexPP2 = D.getNearestPointIndex(a, PP2)[0]-offset
    PP2 = (a[1][0,indexPP2], a[1][1,indexPP2], a[1][2,indexPP2])
    
    indexPP3 = indexPP1 + Nt//2
    PP3 = (a[1][0,indexPP3], a[1][1,indexPP3], a[1][2,indexPP3])
    
    N1 = indexPP2-indexPP1+1; N2 = Nt//2-N1+2
    indexPP4 = indexPP3 + N1-1
    PP4 = (a[1][0,indexPP4], a[1][1,indexPP4], a[1][2,indexPP4])
    
    # Lines
    l1 = D.line(P1, P2, N=N1)
    l2 = D.line(P2, P3, N=N2)
    l3 = D.line(P3, P4, N=N1)
    l4 = D.line(P4, P1, N=N2)

    dist1 = D.getLength(l1)
    p1 = D.line(P1, PP1, N=10)
    dist2 = D.getLength(p1)
    Np = int(dist2/dist1*N1)+1
    
    p1 = D.line(P1, PP1, N=Np)
    p2 = D.line(P2, PP2, N=Np)
    p3 = D.line(P3, PP3, N=Np)
    p4 = D.line(P4, PP4, N=Np)
    
    # subzones
    s1 = T.subzone(a, (indexPP1+1,1,1), (indexPP2+1,1,1))
    s2 = T.subzone(a, (indexPP2+1,1,1), (indexPP3+1,1,1))
    s3 = T.subzone(a, (indexPP3+1,1,1), (indexPP4+1,1,1))
    s4 = T.subzone(a, (indexPP4+1,1,1), (Nt,1,1))
    
    # TFIs
    m = G.TFI([l1,l2,l3,l4])
    m1 = G.TFI([s1, p1, p2, l1])
    m2 = G.TFI([s2, p2, p3, l2])
    m3 = G.TFI([s3, p3, p4, l3])
    m4 = G.TFI([s4, p4, p1, l4])
    return [m,m1,m2,m3,m4]

#==============================================================================
# TFI O (N impair)
# Cherche a optimiser weight et offset.
# Weight: poids pour positionner les 4 points interieurs, si weight est grand
# les points sont plus proche de a.
# Offset: decalage des points sur a par rapport a une projection nearest node
# des points interieurs
#==============================================================================
def TFIO(a, weight=None):
    """Generate a transfinite interpolation mesh for 1 input curve."""
    optWeight = 0; optOffset = 0; optScore = 1.e6
    Nt = a[2]
    if Nt//2 - Nt*0.5 == 0: raise ValueError("TFIO: number of points must be odd.")

    # step: si Nt est tres grand, on diminue le nbre de test pour offseter les points
    step = int(Nt / 2 / 50)+1

    # si weight est donne, on l'utilise sinon plage de recherche de 3 a 10
    if weight is None: weightRange = numpy.arange(3.,10.,1.)
    elif isinstance(weight, float): weightRange = numpy.arange(weight, weight+1)
    elif isinstance(weight, int): weightRange = numpy.arange(weight, weight+1)
    elif len(weight) == 2: weightRange = numpy.arange(weight[0],weight[1],1.)
    elif len(weight) == 3: weightRange = numpy.arange(weight[0],weight[1],weight[2])
    else: weightRange = numpy.arange(3.,10.,1.)

    for j in range(-Nt//4, Nt//4+1, step):
        for i in weightRange:
            try:
                [m,m1,m2,m3,m4] = TFIO__(a, i, j)
                score = quality([m,m1,m2,m3,m4])
                if score < optScore:
                    optWeight = i; optOffset = j; optScore = score
            except: pass
    print('Info: TFIO: resulting weight=%g, offset=%g.'%(optWeight,optOffset))
    print('Info: TFIO: resulting score=%g.'%optScore)
    return TFIO__(a, optWeight, optOffset)

#==============================================================================
# Cree un demi-maillage en O
# a1: straight, a2: round
# Celle qui a le plus de points est prise pour round.
#==============================================================================
def TFIHalfO__(a1, a2, weight, offset=0, tol=1.e-6):
    import Transform as T
    import Geom as D

    Nt1 = a1[2]; Nt2 = a2[2]
    if Nt1 > Nt2: ap = a2; a2 = a1; a1 = ap; Np = Nt2; Nt2 = Nt1; Nt1 = Np
    # le plus long est toujours Nt2

    # Check
    P0 = (a1[1][0,0], a1[1][1,0], a1[1][2,0])
    P00 = (a2[1][0,0], a2[1][1,0], a2[1][2,0])
    if abs(P0[0]-P00[0]) > tol: a2 = T.reorder(a2, (-1,2,3))
    elif abs(P0[1]-P00[1]) > tol: a2 = T.reorder(a2, (-1,2,3))
    elif abs(P0[2]-P00[2]) > tol: a2 = T.reorder(a2, (-1,2,3))
    
    # Round
    w1 = C.array('weight', Nt1, 1, 1); w1 = C.initVars(w1, 'weight', 1.)
    w = C.array('weight', Nt2, 1, 1)
    w = C.initVars(w, 'weight', 1.); w[1][0,0:Nt2//2+1] = weight
    P3 = G.barycenter([a2,a1], [w,w1])
    w = C.initVars(w, 'weight', 1.); w[1][0,Nt2//2:Nt2] = weight
    P4 = G.barycenter([a2,a1], [w,w1])

    # Projection
    b = C.convertArray2Hexa(a2); b = G.close(b)
    PP = C.array('x,y,z', 2,1,1)
    C.setValue(PP, 0, (P3[0], P3[1], P3[2]))
    C.setValue(PP, 1, (P4[0], P4[1], P4[2]))
    PPP = T.projectOrtho(PP, [b])
    PPP = PPP[1]

    PP3 = (PPP[0,0], PPP[1,0], PPP[2,0])
    PP4 = (PPP[0,1], PPP[1,1], PPP[2,1])
    
    indexPP3 = D.getNearestPointIndex(a2, PP3)[0]+offset
    if ((Nt1-Nt2)//2-(Nt1-Nt2)*0.5 == 0 and
        (indexPP3+1)//2-(indexPP3+1)*0.5 != 0): indexPP3 += 1
    if ((Nt1-Nt2)//2-(Nt1-Nt2)*0.5 != 0 and
        (indexPP3+1)//2-(indexPP3+1)*0.5 == 0): indexPP3 += 1
    #if (indexPP3 == 0): indexPP3 = 1
    #elif (indexPP3 == (Nt2-1)//2): indexPP3 += -1
    PP3 = (a2[1][0,indexPP3], a2[1][1,indexPP3], a2[1][2,indexPP3])
    N1 = indexPP3+1
    indexPP4 = Nt2-N1 
    PP4 = (a2[1][0,indexPP4], a2[1][1,indexPP4], a2[1][2,indexPP4])
    N2 = Nt2-2*N1+2

    # Straight
    N3 = (Nt1-N2+2)//2
    
    ind = N3-1
    P1 = (a1[1][0,ind], a1[1][1,ind], a1[1][2,ind])
    P2 = (a1[1][0,Nt1-ind-1], a1[1][1,Nt1-ind-1], a1[1][2,Nt1-ind-1])

    # Lines
    l1 = D.line(P1, P3, N=N1)
    l2 = D.line(P3, P4, N=N2)
    l3 = D.line(P4, P2, N=N1)
    
    p1 = D.line(P3, PP3, N=N3)
    p2 = D.line(P4, PP4, N=N3)

    # subzones
    s1 = T.subzone(a1, (1,1,1), (ind+1,1,1) )
    s2 = T.subzone(a1, (ind+1,1,1), (Nt1-ind,1,1) )
    s3 = T.subzone(a1, (Nt1-ind,1,1), (Nt1,1,1) )
    
    s4 = T.subzone(a2, (1,1,1), (indexPP3+1,1,1) )
    s5 = T.subzone(a2, (indexPP3+1,1,1), (indexPP4+1,1,1) )
    s6 = T.subzone(a2, (indexPP4+1,1,1), (Nt2,1,1) )

    #C.convertArrays2File([l1,l2,l3,p1,p2,s1,s2,s3,s4,s5,s6], 'lines.plt')

    # TFIs
    m = G.TFI([l1,l2,l3,s2])
    m1 = G.TFI([s1, s4, p1, l1])
    m1 = T.reorder(m1, (1,-2,3))
    m2 = G.TFI([s5, p1, p2, l2])
    m3 = G.TFI([s6, p2, s3, l3])

    return [m,m1,m2,m3]

#==============================================================================
# TFI half O (N1 et N2 impairs)
#==============================================================================
def TFIHalfO(a1, a2):
    """Generate a transfinite interpolation mesh for 2 input curves."""
    optWeight = 0; optOffset = 0; optScore = 1.e6
    Nt1 = a1[2]; Nt2 = a2[2]
    if Nt1//2 - Nt1*0.5 == 0 and Nt2//2 - Nt2*0.5 != 0:
        raise ValueError("TFIHalfO: N1 and N2 must be odd.")

    step = int(Nt2 / 4 / 50)+1
    for j in range(-Nt2//8, Nt2//8, step):
        for i in range(2, 10):
            try:
                [m,m1,m2,m3] = TFIHalfO__(a1, a2, i, j)
                score = quality([m,m1,m2,m3])
                if score < optScore:
                    optWeight = i; optScore = score; optOffset = j
            except: pass
    print('Info: TFIHalfO; resulting score=%g'%optScore)
    return TFIHalfO__(a1, a2, optWeight, optOffset)

#==============================================================================
# Cree un seul maillage a partir de 2 courbes
# a1: straight, a2: round
# Il faut N1-N2 pair.
# Celle qui a le plus de points est prise pour round.
#==============================================================================
def TFIMono(a1, a2):
    """Generate a transfinite interpolation mesh for 2 input curves."""
    import Transform as T
    N1 = a1[2]; N2 = a2[2]
    diff = N2-N1
    if diff == 0:
        Np = N1//2
        b1 = T.subzone(a1, (1,1,1), (Np+1,1,1))
        b2 = T.subzone(a1, (Np+1,1,1), (N1,1,1))
        b3 = T.subzone(a2, (1,1,1), (N1-Np,1,1))
        b4 = T.subzone(a2, (N1-Np,1,1), (N2,1,1))
        #C.convertArrays2File([b1,b2,b3,b4], 'test.plt')
        m1 = G.TFI([b1,b2,b3,b4])
        return [m1]

    if diff//2 != diff*0.5: raise ValueError("TFIMono: N1-N2 must be even.")
    if diff < 0: ap = a2; a2 = a1; a1 = ap; diff = -diff; N2 = N1
    Np = (diff+2)//2
    b1 = T.subzone(a2, (1,1,1), (Np,1,1))
    b2 = T.subzone(a2, (N2-Np+1,1,1), (N2,1,1))
    b3 = T.subzone(a2, (Np,1,1), (N2-Np+1,1,1))
    m1 = G.TFI([a1,b1,b3,b2])
    return [m1]

#==============================================================================
# Cree un ensemble de maillages TFI en etoilant les edges et en 
# faisant des TFIs par triangle
#==============================================================================
def TFIStar(edges):
    """Generate a transfinite interpolation mesh for a list of input curves."""
    import Geom as D
    XG = G.barycenter(edges)
    out = []
    for e in edges:
        N = e[2]
        ep = e[1]
        P0 = (ep[0,0], ep[1,0], ep[2,0])
        P1 = (ep[0,-1], ep[1,-1], ep[2,-1])
        l1 = D.line(P0, XG, N=N)
        l2 = D.line(XG, P1, N=N)
        ret = TFITri(e, l1, l2)
        out += ret
    return out

# Cree un ensemble de maillages TFI en etoilant les milieux des edges
# Les edges doivent tous avoir le meme nombre de points impair
def TFIStar2(edges):
    """Generate a transfinite interpolation mesh for a list of input curves."""
    import Geom as D
    import Transform as T
    orderEdges(edges, tol=1.e-6)
    
    XG = G.barycenter(edges) # a optimiser
    out = []
    for c, e in enumerate(edges):
        N1 = e[2]//2+1
        e1 = T.subzone(e, (N1,1,1), (-1,1,1))
        if c == len(edges)-1: en = edges[0]
        else: en = edges[c+1]
        N2 = en[2]//2+1
        e2 = T.subzone(en, (1,1,1), (N2,1,1))
        e1p = e1[1]; e2p = e2[1]
        P0 = (e1p[0,0], e1p[1,0], e1p[2,0])
        P1 = (e2p[0,-1], e2p[1,-1], e2p[2,-1])
        l1 = D.line(P0,XG,N=N2)
        l2 = D.line(XG,P1,N=N1)
        ret = G.TFI([e1,e2,l1,l2])
        #ret = [e1,e2,l1,l2]
        out += [ret]
    
    return out
    