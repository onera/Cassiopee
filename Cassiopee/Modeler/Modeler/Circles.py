"""All models of circles."""
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Cercle avec dents
# IN: R: rayon cercle
# IN: Rd: rayon dents
#==============================================================================
def circle1(R=1., Rd=0.8, Nd=10, fracD=0.5, N=180):
    """Circle"""
    c1 = D.circle((0,0,0), R, N=N)
    c2 = D.circle((0,0,0), Rd, N=N)
    Npts = int(N*1./Nd)
    Npts1 = int((1.-fracD)*Npts)
    out = []; c = 1
    for i in range(Nd):
        bx = min(c+Npts, N)
        a1 = T.subzone(c1, (c,1,1), (c+Npts1,1,1))
        a2 = T.subzone(c2, (c+Npts1,1,1), (bx,1,1))
        a3 = D.line(C.getValue(c1,c+Npts1-1), C.getValue(c2,c+Npts1-1), N=3)
        a4 = D.line(C.getValue(c1,bx-1), C.getValue(c2,bx-1), N=3)
        out += [a1,a2,a3,a4]
        c += Npts
    out = C.convertArray2Hexa(out)
    out = T.join(out)
    out = G.close(out)
    return out

#=============================================================================
# Cercle avec bumps
#=============================================================================
def circle2(R=1., Rd=0.8, Nd=10, fracD=0.5, N=180):
    """Circle"""
    c1 = D.circle((0,0,0), R, N=N)
    c2 = D.circle((0,0,0), Rd, N=N)
    Npts = int(N*1./Nd)
    Npts1 = int((1.-fracD)*Npts)
    out = []; c = 1
    for i in range(Nd):
        bx = min(c+Npts, N)
        a1 = T.subzone(c1, (c,1,1), (c+Npts1,1,1))
        a2 = T.subzone(c2, (c+Npts1,1,1), (bx,1,1))
        a3 = D.line(C.getValue(c1,c+Npts1-1), C.getValue(c2,c+Npts1-1), N=3)
        a4 = D.line(C.getValue(c1,bx-1), C.getValue(c2,bx-1), N=3)
        a2 = T.oneovern(a2, (Npts1-1,1,1))
        p = [a3,a2,a4]
        p = T.join(p)
        s = D.spline(p, N=Npts1)
        out += [a1,s]
        c += Npts
    out = C.convertArray2Hexa(out)
    out = T.join(out)
    out = G.close(out)
    return out
