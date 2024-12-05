"""Offset a surface."""
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Post.PyTree as P
import Dist2Walls.PyTree as DTW

#==============================================================================
# Calcul la distance a la paroi par differents algos
# IN: b: maillage ou l'on calcule la distance
# IN: a: maillage du corps
# IN: loc='centers', 'nodes': localisation du champ de distance
#==============================================================================
def _compDistance__(b, a, loc):
    try:
        DTW._distance2Walls(b, a, type='ortho', loc=loc, signed=1)
        fail = False
    except: fail = True
    if fail:
        try:
            DTW._distance2Walls(b, a, type='ortho', loc=loc, signed=0)
            fail = False
        except: fail = True
    if fail:
        try:
            DTW._distance2Walls(b, a, type='mininterf', loc=loc, signed=0)
            fail = False
        except: raise
    return None

#==============================================================================
def withCart__(a, offset, pointsPerUnitLength, dim=3):
    # Grille cartesienne
    BB = G.bbox(a)
    xmin = BB[0]; ymin = BB[1]; zmin = BB[2]
    xmax = BB[3]; ymax = BB[4]; zmax = BB[5]
    ni = pointsPerUnitLength*(xmax-xmin);
    nj = pointsPerUnitLength*(ymax-ymin)
    nk = pointsPerUnitLength*(zmax-zmin)

    if ni < 2: ni = 2
    if nj < 2: nj = 2
    if nk < 2: nk = 2
    hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1); hk = (zmax-zmin)/(nk-1)
    h = min(hi, hj);
    if dim == 3: h = min(h, hk);
    h = max(h, 1.e-6)
    ni = int((xmax-xmin)/h)+7; nj = int((ymax-ymin)/h)+7
    nk = int((zmax-zmin)/h)+7
    ni += int(2*abs(offset)/h); nj += int(2*abs(offset)/h); nk += int(2*abs(offset)/h)
    xc = (xmax+xmin)*0.5; Lx = ni*h*0.5
    yc = (ymax+ymin)*0.5; Ly = nj*h*0.5
    zc = (zmax+zmin)*0.5; Lz = nk*h*0.5

    if dim == 2:
        nk = 1
        b = G.cart( (xc-Lx, yc-Ly, zmin), (h, h, h), (ni,nj,nk) )
    else:
        b = G.cart( (xc-Lx, yc-Ly, zc-Lz), (h, h, h), (ni,nj,nk) )

    # Calcul la distance a la paroi
    _compDistance__(b, a, loc='nodes')
    # Extraction isoSurf
    iso = P.isoSurfMC([b], 'TurbulentDistance', value=offset)
    C._rmVars(iso,'TurbulentDistance')
    return iso

#==============================================================================
def withOctree__(a, offset, pointsPerUnitLength, dim=3):
    # step
    tol = 1./pointsPerUnitLength

    # octree
    snears = []; sec = 0
    for z in Internal.getZones(a):
        bb = G.bbox(z)
        rx = bb[3]-bb[0]; ry = bb[4]-bb[1]; rz = bb[5]-bb[2]
        snear = min(rx, ry);
        if dim==3: snear = min(snear, rz)
        snear = 0.1*snear
        sec = max(sec, snear)
        snears.append(snear)
    o = G.octree(a, snears, dfar=offset+sec)
    _compDistance__(o, a, loc='nodes')

    # iteration d'adaptation
    nit = 0
    while nit < 10:
        print('iterating: %d...'%nit)

        o = C.node2Center(o, 'TurbulentDistance')
        G._getVolumeMap(o)

        # adapt
        C._initVars(o, '{centers:vol}={centers:vol}**0.33333')
        # was 2.1 factor
        C._initVars(o, '{centers:indicator}=logical_and({centers:vol} > %20.16g , abs({centers:TurbulentDistance}-%20.16g) < 1.*{centers:vol})'%(tol,offset))
        o1 = G.adaptOctree(o, 'centers:indicator')

        #C.convertPyTree2File([o1]+a, 'out%d.cgns'%nit)

        # check convergence
        dim1 = Internal.getZoneDim(o1); dim = Internal.getZoneDim(o)
        if dim1 == dim: break

        #if (nit%2 == 0): o1 = P.extractMesh([o], o1)
        _compDistance__(o1, a, loc='nodes')
        o = o1
        nit += 1

    C._rmVars(o, ['centers:TurbulentDistance','centers:vol','centers:indicator'])
    #C.convertPyTree2File(o, 'out.cgns')

    # Iso surface
    iso = P.isoSurfMC([o], 'TurbulentDistance', value=offset)
    return iso



def offsetSurface(a, offset=1., pointsPerUnitLength=1., algo=0, dim=3):
    """Offset a surface a given distance.
    Usage: offsetSurface(t,offset, pointsPerUnitLength, algo, dim)"""
    if algo==0: iso = withCart__(Internal.getZones(a), offset, pointsPerUnitLength, dim)
    else: iso = withOctree__(Internal.getZones(a), offset, pointsPerUnitLength, dim)
    return iso
