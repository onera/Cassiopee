# - All panels -
import Geom as D
import Transform as T
import Generator as G
import Converter as C
from . import Boxes
from . import Cylinders

#==============================================================================
# Panneau avec un texte
# text: le texte a afficher (\n possible)
# scale: l'echelle
# h: hauteur du poteau
#==============================================================================
def panel(text="Hello", scale=0.1, h=1.):
    t = D.text2D(text, font='vera')
    t = T.rotate(t, (0,0,0), (1,0,0), 90.)
    t = T.scale(t, scale)
    bb = G.bbox(t)
    t = T.translate(t, (0,0,-bb[2]+h))
    bb = G.bbox(t)
    dt = 0.3
    box1 = Boxes.box((bb[0]-2*dt,bb[1]+0.001,bb[2]-2*dt), (bb[3]+2*dt,bb[4]+2*dt,bb[5]+2*dt), chamfer=0.25*dt)
    box1 = C.convertArray2Tetra(box1)
    bb = G.bbox(box1)

    dx = bb[3]-bb[0]; dy = bb[4]-bb[1]; dz = bb[5]-bb[2]
    poteau = Cylinders.cylinder(R1=0.75*dt,R2=0.75*dt,h=h)
    poteau = T.translate(poteau, (bb[0]+dx*0.5,bb[1]+dy*0.5,bb[2]-h))
    o = T.join([t,box1,poteau])
    o = G.close(o)
    return o

#==============================================================================
# Panneau avec une frame carree
# w: width
# h: height
# border: largeur du border
#==============================================================================
def frame(w, h, border):
    P0 = (-w*0.5, -0.5*h, 0.)
    P1 = (w*0.5, -0.5*h, 0.)
    P2 = (w*0.5, +0.5*h, 0.)
    P3 = (-w*0.5, +0.5*h, 0.)
    a = Boxes.box2D(P0, P2, uv=True)

    P4 = (P0[0]-border,P0[1]-border,P0[2]+border)
    P4b = (P0[0]-border,P0[1]-border,P0[2])
    P5 = (P1[0]+border,P1[1]-border,P1[2]+border)
    P5b = (P1[0]+border,P1[1]-border,P1[2])
    P6 = (P2[0]+border,P2[1]+border,P2[2]+border)
    P6b = (P2[0]+border,P2[1]+border,P2[2])
    P7 = (P3[0]-border,P3[1]+border,P3[2]+border)
    P7b = (P3[0]-border,P3[1]+border,P3[2])

    q1 = D.quadrangle(P0,P1,P5,P4)
    q2 = D.quadrangle(P1,P5,P6,P2)
    q3 = D.quadrangle(P3,P2,P6,P7)
    q4 = D.quadrangle(P0,P3,P7,P4)

    q5 = D.quadrangle(P4,P5,P5b,P4b)
    q6 = D.quadrangle(P5,P6,P6b,P5b)
    q7 = D.quadrangle(P6,P7,P7b,P6b)
    q8 = D.quadrangle(P7,P4,P4b,P7b)

    q9 = D.quadrangle(P4b, P5b, P1, P0)
    q10 = D.quadrangle(P5b, P6b, P2, P1)
    q11 = D.quadrangle(P6b, P7b, P3, P2)
    q12 = D.quadrangle(P7b, P4b, P0, P3)

    all = [q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12]
    all = C.initVars(all, '_u_=0')
    all = C.initVars(all, '_v_=0')
    all = C.convertArray2Tetra(all)
    f = T.join(all)

    return [a,f]