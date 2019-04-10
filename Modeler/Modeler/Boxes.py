# - all boxes -
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Boite
# IN: Pmin, Pmax: pts extremes
# Si chamfer>0: chanfrein droit
#==============================================================================
def box(Pmin, Pmax, chamfer=-1.):
    (xmin,ymin,zmin) = Pmin
    (xmax,ymax,zmax) = Pmax
    if chamfer <= 0.: # pas de chanfrein
        P1 = (xmin,ymin,zmin); P2 = (xmax,ymin,zmin)
        P3 = (xmax,ymax,zmin); P4 = (xmin,ymax,zmin)
        P5 = (xmin,ymin,zmax); P6 = (xmax,ymin,zmax)
        P7 = (xmax,ymax,zmax); P8 = (xmin,ymax,zmax)
        Q1 = D.quadrangle(P1,P2,P3,P4)
        Q2 = D.quadrangle(P5,P6,P7,P8)
        Q3 = D.quadrangle(P1,P4,P8,P5)
        Q4 = D.quadrangle(P2,P3,P7,P6)
        Q5 = D.quadrangle(P1,P2,P6,P5)
        Q6 = D.quadrangle(P4,P3,P7,P8)
        a = [Q1,Q2,Q3,Q4,Q5,Q6]
        a = T.join(a)
        a = G.close(a)
        a = T.reorder(a, (1,))
        return a
    else: # chanfrein droit
        deltax = xmax-xmin
        deltay = ymax-ymin
        deltaz = zmax-zmin
        chamfer = min(chamfer, 0.5*deltax, 0.5*deltay, 0.5*deltaz)
        xminp = xmin+chamfer
        yminp = ymin+chamfer
        zminp = zmin+chamfer
        xmaxp = xmax-chamfer
        ymaxp = ymax-chamfer
        zmaxp = zmax-chamfer
        P1 = (xminp,yminp,zmin); P2 = (xmaxp,yminp,zmin)
        P3 = (xmaxp,ymaxp,zmin); P4 = (xminp,ymaxp,zmin)

        P5 = (xminp,yminp,zmax); P6 = (xmaxp,yminp,zmax)
        P7 = (xmaxp,ymaxp,zmax); P8 = (xminp,ymaxp,zmax)

        P9 = (xmin,yminp,zminp); P10 = (xmin,ymaxp,zminp)
        P11 = (xmin,ymaxp,zmaxp); P12 = (xmin,yminp,zmaxp)

        P13 = (xmax,yminp,zminp); P14 = (xmax,ymaxp,zminp)
        P15 = (xmax,ymaxp,zmaxp); P16 = (xmax,yminp,zmaxp)
        
        P17 = (xminp,ymin,zminp); P18 = (xmaxp,ymin,zminp)
        P19 = (xmaxp,ymin,zmaxp); P20 = (xminp,ymin,zmaxp)

        P21 = (xminp,ymax,zminp); P22 = (xmaxp,ymax,zminp)
        P23 = (xmaxp,ymax,zmaxp); P24 = (xminp,ymax,zmaxp)

        Q1 = D.quadrangle(P1,P4,P3,P2)
        Q2 = D.quadrangle(P5,P6,P7,P8)

        Q3 = D.quadrangle(P9,P12,P11,P10)
        Q4 = D.quadrangle(P13,P14,P15,P16)

        Q5 = D.quadrangle(P17,P18,P19,P20)
        Q6 = D.quadrangle(P21,P24,P23,P22)

        Q7 = D.quadrangle(P1,P2,P18,P17)
        Q8 = D.quadrangle(P3,P4,P21,P22)
        Q9 = D.quadrangle(P2,P3,P14,P13)
        Q10 = D.quadrangle(P4,P1,P9,P10)
        
        Q11 = D.quadrangle(P5,P20,P19,P6)
        Q12 = D.quadrangle(P6,P16,P15,P7)
        Q13 = D.quadrangle(P7,P23,P24,P8)
        Q14 = D.quadrangle(P8,P11,P12,P5)

        Q15 = D.quadrangle(P9,P17,P20,P12)
        Q16 = D.quadrangle(P10,P11,P24,P21)
        Q17 = D.quadrangle(P14,P22,P23,P15)
        Q18 = D.quadrangle(P16,P19,P18,P13)

        Q19 = D.quadrangle(P5,P12,P20,P20)
        Q20 = D.quadrangle(P11,P8,P24,P24)
        Q21 = D.quadrangle(P16,P6,P19,P19)
        Q22 = D.quadrangle(P7,P15,P23,P23)

        Q23 = D.quadrangle(P9,P1,P17,P17)
        Q24 = D.quadrangle(P4,P10,P21,P21)
        Q25 = D.quadrangle(P14,P3,P22,P22)
        Q26 = D.quadrangle(P2,P13,P18,P18)

        a = [Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,
             Q18,Q19,Q20,Q21,Q22,Q23,Q24,Q25,Q26]
        a = T.join(a)
        a = G.close(a)
        return a

#==========================================================================
# IN: Pmin, Pmax: pts min et max
# IN: r: % de round (entre 0 et 1)
#==========================================================================
def box2D(Pmin, Pmax, r=0.):
    xmin = Pmin[0]; ymin = Pmin[1]
    xmax = Pmax[0]; ymax = Pmax[1]
    dx = r*(xmax-xmin); dy = r*(ymax-ymin) 
    P0 = (xmin,ymin,0)
    P1 = (xmin,ymin+dy,0)
    P2 = (xmin,ymax-dy,0)
    P3 = (xmin,ymax,0)
    P4 = (xmin+dx,ymax,0)
    P5 = (xmax-dx,ymax,0)
    P6 = (xmax,ymax,0)
    P7 = (xmax,ymax-dy,0)
    P8 = (xmax,ymin+dy,0)
    P9 = (xmax, ymin,0)
    P10 = (xmax-dx,ymin,0)
    P11 = (xmin+dx,ymin,0)
    d1 = D.line(P1, P2, N=4)
    d2 = D.line(P4, P5, N=4)
    d3 = D.line(P7, P8, N=4)
    d4 = D.line(P10, P11, N=4)
    p1 = D.polyline([P11,P0,P1])
    p2 = D.polyline([P2,P3,P4])
    p3 = D.polyline([P5,P6,P7])
    p4 = D.polyline([P8,P9,P10])
    s1 = D.spline(p1, N=10)
    s2 = D.spline(p2, N=10)
    s3 = D.spline(p3, N=10)
    s4 = D.spline(p4, N=10)
    if r > 1.e-12: a = [s1,s2,s3,s4,d1,d2,d3,d4]
    else: a = [d1,d2,d3,d4]
    a = C.convertArray2Tetra(a)
    a = T.join(a)
    a = G.close(a)
    a = G.tetraMesher(a)
    return a

#=================================================================
# Ellipse 2D
#=================================================================
def ellipse2D(Pmin, Pmax):
    xmin = Pmin[0]; ymin = Pmin[1]
    xmax = Pmax[0]; ymax = Pmax[1]
    xc = (xmin+xmax)*0.5
    yc = (ymin+ymax)*0.5
    dx = xmax-xmin
    dy = ymax-ymin
    R = 1.4142*dx*0.5
    a = D.circle((xc,yc,0), R, N=40)
    a = T.contract(a, (xc,yc,0),(1,0,0),(0,0,1),dy*1./dx)
    a = G.tetraMesher(a)
    return a

#=================================================================
# deformNormals ne marche pas...
def ellipseSpikes2D(Pmin, Pmax, r=0.1):
    xmin = Pmin[0]; ymin = Pmin[1]
    xmax = Pmax[0]; ymax = Pmax[1]
    xc = (xmin+xmax)*0.5
    yc = (ymin+ymax)*0.5
    dx = xmax-xmin
    dy = ymax-ymin
    R = 1.4142*dx*0.5
    a = D.circle((xc,yc,0), R, N=20)
    a = C.convertArray2Hexa(a)
    a = G.close(a)
    b = C.initVars(a, 'F', 0.)
    b = C.extractVars(b, ['F'])
    #for i in range(0,20,2): b[1][0,i] = r
    #for i in range(1,20,2): b[1][0,i] = -r
    a = T.deformNormals(a, b)
    a = T.contract(a, (xc,yc,0),(1,0,0),(0,0,1),dy*1./dx)
    a = G.tetraMesher(a)
    return a    
