"""PolyQuad mesh generator. Extension of Generator.
"""
from . import Generator as G
import Geom as D
import math

try: range = xrange
except: pass

__version__ = G.__version__

#=============================================================================
# PolyQuadMesher pour les QUAD-arrays
#=============================================================================
def polyQuadMesher(polyQuad, h, hf, density, next):
    """Generate a multiple mesh for a polyquad.
    Usage:
    polyQuadMesher(polyQuad, h, hf, density, next)"""
    import Converter as C

    polyQuad = G.close(polyQuad)

    addFactor = 0.2
    if len(polyQuad) != 4:
        raise TypeError("polyQuadMesher: requires a QUAD array.")
    else:
        if polyQuad[3] != 'QUAD':
            raise TypeError("polyQuadMesher: requires a QUAD array.")

    f = polyQuad[1]; c = polyQuad[2]; ne = c.shape[1]

    deuxPiSur3 = 2.*math.pi/3.

    # Calcul des longueurs minimum et maximum des arretes
    lmin = 1.e6
    lmax = 0.
    for i in range(ne):
        ind1 = c[0,i]-1; ind2 = c[1,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        lmin = min(lmin, l)
        lmax = max(lmax, l)
        ind1 = c[1,i]-1; ind2 = c[2,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        lmin = min(lmin, l)
        lmax = max(lmax, l)
        ind1 = c[2,i]-1; ind2 = c[3,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        lmin = min(lmin, l)
        lmax = max(lmax, l)
        ind1 = c[3,i]-1; ind2 = c[0,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        l = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
        lmin = min(lmin, l)
        lmax = max(lmax, l)

    # Detection de la hauteur maximum admissible
    if h > 0.9*lmin:
        h = 0.9*lmin
        print("Warning: height changed to", h,"...")
        print("...because length of a line segment is", lmin)

    # Detection de la densite minimum
    nk = int(h*density)+1
    if nk < 4:
        density = 4./h
        print("Warning: density changed to", density, "to have 4 points in height.")
    n = int(lmax*density)+1
    if n < 4:
        density = 4./lmax
        print("Warning: density changed to", density, "to have 4 points per segment",i)

    # Calcul automatique de l'extension
#    extension = int(h*density)

    # Calculs prealables
    nk = int(h*density)+1
    distribk = G.cart((0,0,0), (1./nk,1,1), (nk+1,1,1))
    add = max(nk * h / (20 * hf), 1)
    add = min(int(add), 3*nk)
    add = int(addFactor*add)
    distribk = G.enforcePlusX(distribk, hf/h, nk-1, add)
    nk = distribk[2]-1
    delta = C.array('d',nk,1,1)
    for i in range(nk):
        delta[1][0,i] = h*( distribk[1][0,i+1] - distribk[1][0,i])
    mesh = []; walls = []

    # Generation des maillages
    for i in range(ne):
        ind1 = c[0,i]-1; ind2 = c[1,i]-1 ; ind3 = c[2,i]-1; ind4 = c[3,i]-1
        x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
        x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
        x3 = f[0,ind3]; y3 = f[1,ind3]; z3 = f[2,ind3]
        x4 = f[0,ind4]; y4 = f[1,ind4]; z4 = f[2,ind4]

        # ext = 1 ; extension chimere
        # ext = 0 ; TFI paroi
        # ext = -1; TFI MD + TTM
        iQ1 = findNeighbourIndex(polyQuad,i,ind1,ind2)
        iQ2 = findNeighbourIndex(polyQuad,i,ind2,ind3)
        iQ3 = findNeighbourIndex(polyQuad,i,ind3,ind4)
        iQ4 = findNeighbourIndex(polyQuad,i,ind4,ind1)
        [nx,ny,nz] = normalVector(polyQuad,i)
        [n1x,n1y,n1z] = normalVector(polyQuad,iQ1)
        [n2x,n2y,n2z] = normalVector(polyQuad,iQ2)
        [n3x,n3y,n3z] = normalVector(polyQuad,iQ3)
        [n4x,n4y,n4z] = normalVector(polyQuad,iQ4)
        [t1x,t1y,t1z] = tangentVector(polyQuad,i,ind1,ind2)
        [t2x,t2y,t2z] = tangentVector(polyQuad,i,ind2,ind3)
        [t3x,t3y,t3z] = tangentVector(polyQuad,i,ind3,ind4)
        [t4x,t4y,t4z] = tangentVector(polyQuad,i,ind4,ind1)

        a1 = -(nx*n1x + ny*n1y + nz*n1z)
        a1 = min(a1,1.); a1 = max(a1,-1.)
        a2 = -(nx*n2x + ny*n2y + nz*n2z)
        a2 = min(a2,1.); a2 = max(a2,-1.)
        a3 = -(nx*n3x + ny*n3y + nz*n3z)
        a3 = min(a3,1.); a3 = max(a3,-1.)
        a4 = -(nx*n4x + ny*n4y + nz*n4z)
        a4 = min(a4,1.); a4 = max(a4,-1.)
        if t1x*n1x + t1y*n1y + t1z*n1z < -1.e-10: # extension
            ext1 = 1
        elif math.acos(a1) < deuxPiSur3: # TFI 2 parois
            ext1 = 0
        else: # TFI MD + TTM
            ext1 = -1
        if t2x*n2x + t2y*n2y + t2z*n2z < -1.e-10: # extension
            ext2 = 1
        elif math.acos(a2) < deuxPiSur3: # TFI 2 parois
            ext2 = 0
        else: # TFI MD + TTM
            ext2 = -1
        if t3x*n3x + t3y*n3y + t3z*n3z < -1.e-10: # extension
            ext3 = 1
        elif math.acos(a3) < deuxPiSur3: # TFI 2 parois
            ext3 = 0
        else: # TFI MD + TTM
            ext3 = -1
        if t4x*n4x + t4y*n4y + t4z*n4z < -1.e-10: # extension
            ext4 = 1
        elif math.acos(a4) < deuxPiSur3: # TFI 2 parois
            ext4 = 0
        else: # TFI MD + TTM
            ext4 = -1

        lext1 = max(0,ext1)
        lext2 = max(0,ext2)
        lext3 = max(0,ext3)
        lext4 = max(0,ext4)
        n = int(lmax*density)+1
        ni = n + (lext2+lext4)*next
        nj = n + (lext1+lext3)*next
        nextsn1 = next/(n-1.)

        # p1
        q0x = x4 + lext4*(x4-x3)*nextsn1
        q0y = y4 + lext4*(y4-y3)*nextsn1
        q0z = z4 + lext4*(z4-z3)*nextsn1
        q1x = x1 + lext4*(x1-x2)*nextsn1
        q1y = y1 + lext4*(y1-y2)*nextsn1
        q1z = z1 + lext4*(z1-z2)*nextsn1
        r0x = x2 + lext1*(x2-x3)*nextsn1
        r0y = y2 + lext1*(y2-y3)*nextsn1
        r0z = z2 + lext1*(z2-z3)*nextsn1
        r1x = x1 + lext1*(x1-x4)*nextsn1
        r1y = y1 + lext1*(y1-y4)*nextsn1
        r1z = z1 + lext1*(z1-z4)*nextsn1

        ux = q1x - q0x
        uy = q1y - q0y
        uz = q1z - q0z
        vx = (r1y - r0y)*nz - (r1z - r0z)*ny
        vy = (r1z - r0z)*nx - (r1x - r0x)*nz
        vz = (r1x - r0x)*ny - (r1y - r0y)*nx
        wx = q0x - r0x
        wy = q0y - r0y
        wz = q0z - r0z
        s = - (vx*wx + vy*wy + vz*wz) / (vx*ux + vy*uy + vz*uz)

        p1x = q0x + s*(q1x-q0x)
        p1y = q0y + s*(q1y-q0y)
        p1z = q0z + s*(q1z-q0z)

        # p2
        q0x = x1 + lext1*(x1-x4)*nextsn1
        q0y = y1 + lext1*(y1-y4)*nextsn1
        q0z = z1 + lext1*(z1-z4)*nextsn1
        q1x = x2 + lext1*(x2-x3)*nextsn1
        q1y = y2 + lext1*(y2-y3)*nextsn1
        q1z = z2 + lext1*(z2-z3)*nextsn1
        r0x = x3 + lext2*(x3-x4)*nextsn1
        r0y = y3 + lext2*(y3-y4)*nextsn1
        r0z = z3 + lext2*(z3-z4)*nextsn1
        r1x = x2 + lext2*(x2-x1)*nextsn1
        r1y = y2 + lext2*(y2-y1)*nextsn1
        r1z = z2 + lext2*(z2-z1)*nextsn1

        ux = q1x - q0x
        uy = q1y - q0y
        uz = q1z - q0z
        vx = (r1y - r0y)*nz - (r1z - r0z)*ny
        vy = (r1z - r0z)*nx - (r1x - r0x)*nz
        vz = (r1x - r0x)*ny - (r1y - r0y)*nx
        wx = q0x - r0x
        wy = q0y - r0y
        wz = q0z - r0z
        s = - (vx*wx + vy*wy + vz*wz) / (vx*ux + vy*uy + vz*uz)

        p2x = q0x + s*(q1x-q0x)
        p2y = q0y + s*(q1y-q0y)
        p2z = q0z + s*(q1z-q0z)

        # p3
        q0x = x2 + lext2*(x2-x1)*nextsn1
        q0y = y2 + lext2*(y2-y1)*nextsn1
        q0z = z2 + lext2*(z2-z1)*nextsn1
        q1x = x3 + lext2*(x3-x4)*nextsn1
        q1y = y3 + lext2*(y3-y4)*nextsn1
        q1z = z3 + lext2*(z3-z4)*nextsn1
        r0x = x4 + lext3*(x4-x1)*nextsn1
        r0y = y4 + lext3*(y4-y1)*nextsn1
        r0z = z4 + lext3*(z4-z1)*nextsn1
        r1x = x3 + lext3*(x3-x2)*nextsn1
        r1y = y3 + lext3*(y3-y2)*nextsn1
        r1z = z3 + lext3*(z3-z2)*nextsn1

        ux = q1x - q0x
        uy = q1y - q0y
        uz = q1z - q0z
        vx = (r1y - r0y)*nz - (r1z - r0z)*ny
        vy = (r1z - r0z)*nx - (r1x - r0x)*nz
        vz = (r1x - r0x)*ny - (r1y - r0y)*nx
        wx = q0x - r0x
        wy = q0y - r0y
        wz = q0z - r0z
        s = - (vx*wx + vy*wy + vz*wz) / (vx*ux + vy*uy + vz*uz)

        p3x = q0x + s*(q1x-q0x)
        p3y = q0y + s*(q1y-q0y)
        p3z = q0z + s*(q1z-q0z)

        # p4
        q0x = x3 + lext3*(x3-x2)*nextsn1
        q0y = y3 + lext3*(y3-y2)*nextsn1
        q0z = z3 + lext3*(z3-z2)*nextsn1
        q1x = x4 + lext3*(x4-x1)*nextsn1
        q1y = y4 + lext3*(y4-y1)*nextsn1
        q1z = z4 + lext3*(z4-z1)*nextsn1
        r0x = x1 + lext4*(x1-x2)*nextsn1
        r0y = y1 + lext4*(y1-y2)*nextsn1
        r0z = z1 + lext4*(z1-z2)*nextsn1
        r1x = x4 + lext4*(x4-x3)*nextsn1
        r1y = y4 + lext4*(y4-y3)*nextsn1
        r1z = z4 + lext4*(z4-z3)*nextsn1

        ux = q1x - q0x
        uy = q1y - q0y
        uz = q1z - q0z
        vx = (r1y - r0y)*nz - (r1z - r0z)*ny
        vy = (r1z - r0z)*nx - (r1x - r0x)*nz
        vz = (r1x - r0x)*ny - (r1y - r0y)*nx
        wx = q0x - r0x
        wy = q0y - r0y
        wz = q0z - r0z
        s = - (vx*wx + vy*wy + vz*wz) / (vx*ux + vy*uy + vz*uz)

        p4x = q0x + s*(q1x-q0x)
        p4y = q0y + s*(q1y-q0y)
        p4z = q0z + s*(q1z-q0z)

        if ext1 == 1:
            n1x = ny*(p2z - p1z) - nz*(p2y - p1y)
            n1y = nz*(p2x - p1x) - nx*(p2z - p1z)
            n1z = nx*(p2y - p1y) - ny*(p2x - p1x)
        elif ext1 == -1:
            rx = (n1y+ny)*(z2 - z1) - (n1z+nz)*(y2 - y1)
            ry = (n1z+nz)*(x2 - x1) - (n1x+nx)*(z2 - z1)
            rz = (n1x+nx)*(y2 - y1) - (n1y+ny)*(x2 - x1)
            norme = math.sqrt(rx*rx + ry*ry + rz*rz)
            n1x = rx/norme
            n1y = ry/norme
            n1z = rz/norme

        if ext2 == 1:
            n2x = ny*(p3z - p2z) - nz*(p3y - p2y)
            n2y = nz*(p3x - p2x) - nx*(p3z - p2z)
            n2z = nx*(p3y - p2y) - ny*(p3x - p2x)
        elif ext2 == -1:
            rx = (n2y+ny)*(z3 - z2) - (n2z+nz)*(y3 - y2)
            ry = (n2z+nz)*(x3 - x2) - (n2x+nx)*(z3 - z2)
            rz = (n2x+nx)*(y3 - y2) - (n2y+ny)*(x3 - x2)
            norme = math.sqrt(rx*rx + ry*ry + rz*rz)
            n2x = rx/norme
            n2y = ry/norme
            n2z = rz/norme

        if ext3 == 1:
            n3x = ny*(p4z - p3z) - nz*(p4y - p3y)
            n3y = nz*(p4x - p3x) - nx*(p4z - p3z)
            n3z = nx*(p4y - p3y) - ny*(p4x - p3x)
        elif ext3 == -1:
            rx = (n3y+ny)*(z4 - z3) - (n3z+nz)*(y4 - y3)
            ry = (n3z+nz)*(x4 - x3) - (n3x+nx)*(z4 - z3)
            rz = (n3x+nx)*(y4 - y3) - (n3y+ny)*(x4 - x3)
            norme = math.sqrt(rx*rx + ry*ry + rz*rz)
            n3x = rx/norme
            n3y = ry/norme
            n3z = rz/norme

        if ext4 == 1:
            n4x = ny*(p1z - p4z) - nz*(p1y - p4y)
            n4y = nz*(p1x - p4x) - nx*(p1z - p4z)
            n4z = nx*(p1y - p4y) - ny*(p1x - p4x)
        elif ext4 == -1:
            rx = (n4y+ny)*(z1 - z4) - (n4z+nz)*(y1 - y4)
            ry = (n4z+nz)*(x1 - x4) - (n4x+nx)*(z1 - z4)
            rz = (n4x+nx)*(y1 - y4) - (n4y+ny)*(x1 - x4)
            norme = math.sqrt(rx*rx + ry*ry + rz*rz)
            n4x = rx/norme
            n4y = ry/norme
            n4z = rz/norme

        dh = nx*(x1 + h*nx) + ny*(y1 + h*ny) + nz*(z1 + h*nz)
        d1 = n1x*p1x + n1y*p1y + n1z*p1z
        d2 = n2x*p2x + n2y*p2y + n2z*p2z
        d3 = n3x*p3x + n3y*p3y + n3z*p3z
        d4 = n4x*p4x + n4y*p4y + n4z*p4z

        n41 = nx*n4y*n1z - nx*n4z*n1y + ny*n4z*n1x - ny*n4x*n1z + nz*n4x*n1y - nz*n4y*n1x
        p5x = (dh*(n4y*n1z-n4z*n1y) + d4*(n1y*nz-n1z*ny) + d1*(ny*n4z-nz*n4y))/n41
        p5y = (dh*(n4z*n1x-n4x*n1z) + d4*(n1z*nx-n1x*nz) + d1*(nz*n4x-nx*n4z))/n41
        p5z = (dh*(n4x*n1y-n4y*n1x) + d4*(n1x*ny-n1y*nx) + d1*(nx*n4y-ny*n4x))/n41

        n12 = nx*n1y*n2z - nx*n1z*n2y + ny*n1z*n2x - ny*n1x*n2z + nz*n1x*n2y - nz*n1y*n2x
        p6x = (dh*(n1y*n2z-n1z*n2y) + d1*(n2y*nz-n2z*ny) + d2*(ny*n1z-nz*n1y))/n12
        p6y = (dh*(n1z*n2x-n1x*n2z) + d1*(n2z*nx-n2x*nz) + d2*(nz*n1x-nx*n1z))/n12
        p6z = (dh*(n1x*n2y-n1y*n2x) + d1*(n2x*ny-n2y*nx) + d2*(nx*n1y-ny*n1x))/n12

        n23 = nx*n2y*n3z - nx*n2z*n3y + ny*n2z*n3x - ny*n2x*n3z + nz*n2x*n3y - nz*n2y*n3x
        p7x = (dh*(n2y*n3z-n2z*n3y) + d2*(n3y*nz-n3z*ny) + d3*(ny*n2z-nz*n2y))/n23
        p7y = (dh*(n2z*n3x-n2x*n3z) + d2*(n3z*nx-n3x*nz) + d3*(nz*n2x-nx*n2z))/n23
        p7z = (dh*(n2x*n3y-n2y*n3x) + d2*(n3x*ny-n3y*nx) + d3*(nx*n2y-ny*n2x))/n23

        n34 = nx*n3y*n4z - nx*n3z*n4y + ny*n3z*n4x - ny*n3x*n4z + nz*n3x*n4y - nz*n3y*n4x
        p8x = (dh*(n3y*n4z-n3z*n4y) + d3*(n4y*nz-n4z*ny) + d4*(ny*n3z-nz*n3y))/n34
        p8y = (dh*(n3z*n4x-n3x*n4z) + d3*(n4z*nx-n4x*nz) + d4*(nz*n3x-nx*n3z))/n34
        p8z = (dh*(n3x*n4y-n3y*n4x) + d3*(n4x*ny-n4y*nx) + d4*(nx*n3y-ny*n3x))/n34

        l1 = math.sqrt( (p1x-p2x)*(p1x-p2x) + (p1y-p2y)*(p1y-p2y) + (p1z-p2z)*(p1z-p2z) )
        l2 = math.sqrt( (p2x-p3x)*(p2x-p3x) + (p2y-p3y)*(p2y-p3y) + (p2z-p3z)*(p2z-p3z) )
        l3 = math.sqrt( (p3x-p4x)*(p3x-p4x) + (p3y-p4y)*(p3y-p4y) + (p3z-p4z)*(p3z-p4z) )
        l4 = math.sqrt( (p4x-p1x)*(p4x-p1x) + (p4y-p1y)*(p4y-p1y) + (p4z-p1z)*(p4z-p1z) )

        distribi1 = G.cart((0,0,0), (1./(ni-1),1,1), (ni,1,1))
        distribi2 = G.cart((0,0,0), (1./(ni-1),1,1), (ni,1,1))
        if ext4 == 0:
            distribi1 = G.enforcePlusX(distribi1, hf/max(l1,l3), nk-1, add)
            distribi2 = G.enforceMoinsX(distribi2, hf/max(l1,l3), nk-1, add)
        if ext2 == 0:
            distribi1 = G.enforceMoinsX(distribi1, hf/max(l1,l3), nk-1, add)
            distribi2 = G.enforcePlusX(distribi2, hf/max(l1,l3), nk-1, add)
        distribj = G.cart((0,0,0), (1./(nj-1),1,1), (nj,1,1))
        if ext1 == 0:
            distribj = G.enforcePlusX(distribj, hf/max(l2,l4), nk-1, add)
        if ext3 == 0:
            distribj = G.enforceMoinsX(distribj, hf/max(l2,l4), nk-1, add)

        Q0 = meshQuad((p1x,p1y,p1z), (p2x,p2y,p2z), (p3x,p3y,p3z), (p4x,p4y,p4z), distribi1, distribj)
        Q1 = meshQuad((p2x,p2y,p2z), (p1x,p1y,p1z), (p5x,p5y,p5z), (p6x,p6y,p6z), distribi2, distribk)
        Q2 = meshQuad((p2x,p2y,p2z), (p3x,p3y,p3z), (p7x,p7y,p7z), (p6x,p6y,p6z), distribj, distribk)
        Q3 = meshQuad((p3x,p3y,p3z), (p4x,p4y,p4z), (p8x,p8y,p8z), (p7x,p7y,p7z), distribi2, distribk)
        Q4 = meshQuad((p1x,p1y,p1z), (p4x,p4y,p4z), (p8x,p8y,p8z), (p5x,p5y,p5z), distribj, distribk)
        Qh = meshQuad((p5x,p5y,p5z), (p6x,p6y,p6z), (p7x,p7y,p7z), (p8x,p8y,p8z), distribi1, distribj)

        m = G.TFI([Q4, Q2, Q1, Q3, Q0, Qh])

        mesh.append(m)
        # Walls
        rangesw = []
        if ext4 != 1:
            i1 = 1
        else:
            i1 = next+1
        if ext2 != 1:
            i2 = m[2]
        else:
            i2 = m[2]-next
        if ext1 != 1:
            j1 = 1
        else:
            j1 = next+1
        if ext3 != 1:
            j2 = m[3]
        else:
            j2 = m[3]-next
        wrange = [i1,i2, j1, j2, 1, 1]
        rangesw.append(wrange)

        if ext1 == 0:
            if ext4 == 1:
                if ext2 == 1:
                    wrange = [next+1, m[2]-next, 1, 1, 1, m[4]]
                else:
                    wrange = [next+1, m[2], 1, 1, 1, m[4]]
            else:
                if ext2 == 1:
                    wrange = [1, m[2]-next, 1, 1, 1, m[4]]
                else:
                    wrange = [1, m[2], 1, 1, 1, m[4]]
            rangesw.append(wrange)

        if ext2 == 0:
            if ext1 == 1:
                if ext3 == 1:
                    wrange = [m[2], m[2], next+1, m[3]-next, 1, m[4]]
                else:
                    wrange = [m[2], m[2], next+1, m[3], 1, m[4]]
            else:
                if ext3 == 1:
                    wrange = [m[2], m[2], 1, m[3]-next, 1, m[4]]
                else:
                    wrange = [m[2], m[2], 1, m[3], 1, m[4]]
            rangesw.append(wrange)

        if ext3 == 0:
            if ext4 == 1:
                if ext2 == 1:
                    wrange = [next+1, m[2]-next, m[3], m[3], 1, m[4]]
                else:
                    wrange = [next+1, m[2], m[3], m[3], 1, m[4]]
            else:
                if ext2 == 1:
                    wrange = [1, m[2]-next, m[3], m[3], 1, m[4]]
                else:
                    wrange = [1, m[2], m[3], m[3], 1, m[4]]
            rangesw.append(wrange)

        if ext4 == 0:
            if ext1 == 1:
                if ext3 == 1:
                    wrange = [1, 1, next+1, m[3]-next, 1, m[4]]
                else:
                    wrange = [1, 1, next+1, m[3], 1, m[4]]
            else:
                if ext3 == 1:
                    wrange = [1, 1, 1, m[3]-next, 1, m[4]]
                else:
                    wrange = [1, 1, 1, m[3], 1, m[4]]
            rangesw.append(wrange)

        walls.append(rangesw)

    return [mesh, walls, h, density]

#=============================================================================
# Mesh a quad
#def meshQuad(P1, P2, P3, P4, ni, nj):
#    l1 = D.line(P1,P4,nj)
#    l2 = D.line(P2,P3,nj)
#    l3 = D.line(P1,P2,ni)
#    l4 = D.line(P4,P3,ni)
#    m = G.TFI([l1,l2,l3,l4])
#    return m
#
#=============================================================================
# Mesh a quad
def meshQuad(P1, P2, P3, P4, distrib1, distrib2):
    l1 = D.line(P1,P4)
    l2 = D.line(P2,P3)
    l1 = G.map(l1, distrib2)
    l2 = G.map(l2, distrib2)
    l3 = D.line(P1,P2)
    l4 = D.line(P4,P3)
    l3 = G.map(l3, distrib1)
    l4 = G.map(l4, distrib1)
    m = G.TFI([l1,l2,l3,l4])
    return m

#=============================================================================
# Return the index of neighbour of element i on a QUAD
def findNeighbourIndex(polyQuad,i,iP1,iP2):
    c = polyQuad[2]
    ne = c.shape[1]
    for j in range(ne):
        ind1 = c[0,j]-1; ind2 = c[1,j]-1 ; ind3 = c[2,j]-1; ind4 = c[3,j]-1
        if ind1 == iP1 and ind4 == iP2:
            return j
#            return [j,ind2,ind3]
        if ind2 == iP1 and ind1 == iP2:
            return j
#            return [j,ind3,ind4]
        if ind3 == iP1 and ind2 == iP2:
            return j
#            return [j,ind4,ind1]
        if ind4 == iP1 and ind3 == iP2:
            return j
#            return [j,ind1,ind2]

#=============================================================================
# Return the normal vector to an element i
def normalVector(polyQuad,i):
    c = polyQuad[2]
    f = polyQuad[1]
    ind1 = c[0,i]-1; ind2 = c[1,i]-1 ; ind3 = c[2,i]-1; ind4 = c[3,i]-1
    x1 = f[0,ind1]; y1 = f[1,ind1]; z1 = f[2,ind1]
    x2 = f[0,ind2]; y2 = f[1,ind2]; z2 = f[2,ind2]
    x3 = f[0,ind3]; y3 = f[1,ind3]; z3 = f[2,ind3]
    x4 = f[0,ind4]; y4 = f[1,ind4]; z4 = f[2,ind4]

    nx = (y3 - y1)*(z4 - z2) - (z3 - z1)*(y4 - y2)
    ny = (z3 - z1)*(x4 - x2) - (x3 - x1)*(z4 - z2)
    nz = (x3 - x1)*(y4 - y2) - (y3 - y1)*(x4 - x2)
    norme = math.sqrt(nx*nx + ny*ny + nz*nz)

    if norme == 0.:
        print(i)
        print(ind1, x1, y1, z1)
        print(ind2, x2, y2, z2)
        print(ind3, x3, y3, z3)
        print(ind4, x4, y4, z4)
        raise TypeError("Division par 0! (1)")
    else:
        normi = 1./norme
        nx = nx*normi
        ny = ny*normi
        nz = nz*normi

    return [nx,ny,nz]

#=============================================================================
# Return the tangent vector to an element i directed from the vector P1P2
def tangentVector(polyQuad,i,iP1,iP2):
    f = polyQuad[1]
    x1 = f[0,iP1]; y1 = f[1,iP1]; z1 = f[2,iP1]
    x2 = f[0,iP2]; y2 = f[1,iP2]; z2 = f[2,iP2]

    [nx,ny,nz] = normalVector(polyQuad,i)

    tx = ny*(z2 - z1) - nz*(y2 - y1)
    ty = nz*(x2 - x1) - nx*(z2 - z1)
    tz = nx*(y2 - y1) - ny*(x2 - x1)
    norme = math.sqrt(tx*tx + ty*ty + tz*tz)

    if norme == 0.:
        raise TypeError("Division par 0! (2)")
    else:
        normi = 1./norme
        tx = tx*normi
        ty = ty*normi
        tz = tz*normi

    return [tx,ty,tz]
