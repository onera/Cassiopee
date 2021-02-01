"""All models of tunnels."""
# Fait a partir d'une ligne "x" et de profils empiles
import Generator
import Geom
import Transform
import Converter
import KCore.Vector as Vector

#==================================================================
# IN: les profiles doivent etre dans le plan x,y et centres en 0
# IN: line est une ligne structuree
#==================================================================
def drive(line, profiles):
    n = line[2]
    for i in range(n):
        P = Converter.getValue(line, i)
        if i == n-1:
            PP = Converter.getValue(line, i-1)
            v1 = [P[0]-PP[0],P[1]-PP[1],P[2]-PP[2]]
        else:
            PP = Converter.getValue(line, i+1)
            v1 = [PP[0]-P[0],PP[1]-P[1],PP[2]-P[2]]
        v1 = Vector.normalize(v1)
        v2 = [1,0,0]
        # Partie de v2 ortho a v1
        s = Vector.dot(v1,v2)
        q = Vector.mul(s, v1)
        v2 = Vector.sub(v2, q)
        v2 = Vector.normalize(v2)
        v3 = Vector.cross(v1,v2)
        c = profiles[i]
        c = Transform.rotate(c, (0,0,0), ((1,0,0), (0,1,0), (0,0,1)), (v2,v3,v1))
        profiles[i] = Transform.translate(c, (P[0],P[1],P[2]))
    b = Generator.stack(profiles)
    return b

# Tunnel avec profils circulaires
def tunnel1(line):
    n = line[2]
    profiles = []
    for i in range(n):
        c = Geom.circle((0,0,0),1.,N=30)
        profiles.append(c)
    tunnel = drive(line, profiles)
    return tunnel
