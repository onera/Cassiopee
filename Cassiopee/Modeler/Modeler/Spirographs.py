# - spirograph shapes -
import Geom as D

#=============================================================================
# Courbes spirographes
# IN: 0<=k<=1 ; 0<=l<=1
# IN: Ns: nbre de courbes
# IN: N: nbre de pts sur une courbe
#=============================================================================
def spirograph(k, l, Ns=1, N=100):
    out = []; s = 0.
    for i in range(Ns):
        #print('{x}=(1.-%g)*cos(2*pi*({t}+%g))+%g*%g*cos( (1-%g)/%g*2*pi*({t}+%g); {y}=(1.-%g)*sin(2*pi*({t}+%g))-%g*%g*sin( (1.-%g)/%g*2*pi*({t}+%g)); {z}=0.'%(k,s,k,l,k,k,s,k,s,k,l,k,k,s))
        a = D.curve('{x}=(1.-%g)*cos(2*pi*({t}+%g))+%g*%g*cos( (1-%g)/%g*2*pi*({t}+%g)); {y}=(1.-%g)*sin(2*pi*({t}+%g))-%g*%g*sin( (1.-%g)/%g*2*pi*({t}+%g)); {z}=0.'%(k,s,k,l,k,k,s,k,s,k,l,k,k,s), N)
        #if s*(1.-k)/k > 1.: go = False
        s += 1.
        out.append(a)
    return out
