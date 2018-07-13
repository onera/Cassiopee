"""PolyTri mesh generator. Extension of Generator.
"""
import Generator as G
__version__ = G.__version__

#=============================================================================
# PolyTriMesher pour les TRI-arrays
#=============================================================================
def polyTriMesher(polyTri, h, hf, density, next):
    """Generate a multiple mesh for a polytri.
    Usage:
    polyTriMesher(polyTri, h, hf, density, next)"""
    import PolyQuad
    
    polyQuad = polyTri2PolyQuad(polyTri)
    return PolyQuad.polyQuadMesher(polyQuad, h, hf, density, next)

#=============================================================================
# PolyTri2PolyQuad pour les TRI-arrays
#=============================================================================
def polyTri2PolyQuad(polyTri):
    """Transform a polytri into a polyquad.
    Usage:
    polyTri2PolyQuad(polyTri)"""
    import numpy
    
    polyTri = G.close(polyTri)
    
    if (len(polyTri) != 4):
        raise TypeError("polyTri2PolyQuad: requires a TRI array.")
    else:
        if (polyTri[3] != 'TRI'):
            raise TypeError("polyTri2PolyQuad: requires a TRI array.")
    
    f = polyTri[1]; c = polyTri[2]; ne = c.shape[1]; n = f.shape[1]
    fq = numpy.zeros((3,n+4*ne),dtype=numpy.float64)
    for i in xrange(n):
        fq[0,i] = f[0,i]; fq[1,i] = f[1,i]; fq[2,i] = f[2,i] 
    cq = numpy.zeros((4,3*ne),dtype=numpy.int32)
    
    for i in xrange(ne):
        ind1 = c[0,i]-1; ind2 = c[1,i]-1; ind3 = c[2,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        x3 = f[0,ind3]; y3 = f[1,ind3]; z3 = f[2,ind3]
        
        x4 = (x1+x2+x3)/3
        y4 = (y1+y2+y3)/3
        z4 = (z1+z2+z3)/3

        x5 = (x1+x2)/2
        y5 = (y1+y2)/2
        z5 = (z1+z2)/2

        x6 = (x2+x3)/2
        y6 = (y2+y3)/2
        z6 = (z2+z3)/2

        x7 = (x1+x3)/2
        y7 = (y1+y3)/2
        z7 = (z1+z3)/2
        
        fq[0,n+4*i] = x4; fq[1,n+4*i] = y4; fq[2,n+4*i] = z4
        fq[0,n+4*i+1] = x5; fq[1,n+4*i+1] = y5; fq[2,n+4*i+1] = z5
        fq[0,n+4*i+2] = x6; fq[1,n+4*i+2] = y6; fq[2,n+4*i+2] = z6
        fq[0,n+4*i+3] = x7; fq[1,n+4*i+3] = y7; fq[2,n+4*i+3] = z7

        cq[0,3*i] = ind1+1
        cq[1,3*i] = n+4*i+2
        cq[2,3*i] = n+4*i+1
        cq[3,3*i] = n+4*i+4
        cq[0,3*i+1] = ind2+1
        cq[1,3*i+1] = n+4*i+3
        cq[2,3*i+1] = n+4*i+1
        cq[3,3*i+1] = n+4*i+2
        cq[0,3*i+2] = ind3+1
        cq[1,3*i+2] = n+4*i+4
        cq[2,3*i+2] = n+4*i+1
        cq[3,3*i+2] = n+4*i+3
        
    polyQuad = [polyTri[0], fq, cq, 'QUAD']
    polyQuad = G.close(polyQuad)
    return polyQuad
