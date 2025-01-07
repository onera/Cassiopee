"""Simple Vector module."""
import math

# Multiplie par un scalaire
def mul(scal, v):
    """Multiply a vector by a scalar."""
    return [scal*v[0], scal*v[1], scal*v[2]]

# Ajoute 2 vecteurs (v1+v2)
def add(v1, v2):
    """Add two vectors."""
    return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]

# Soustrait 2 vecteurs (v1-v2)
def sub(v1, v2):
    """Substract two vectors."""
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

# Produit scalaire de 2 vecteurs (v1.v2)
def dot(v1, v2):
    """Scalar product of two vectors."""
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

# Produit vectoriel de 2 vecteurs (v1xv2)
def cross(v1, v2):
    """Vectorial product of two vectors."""
    return [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v2[2]*v1[0], v1[0]*v2[1]-v1[1]*v2[0]]

# Norme d'un vecteur
def norm(v):
    """Vector norm."""
    n = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
    n = math.sqrt(n)
    return n

# Norme au carre d'un vecteur
def norm2(v):
    """Square norm of vector."""
    n = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
    return n

# Normalise le vecteur
def normalize(v):
    """Normalize a vector."""
    n = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
    n = math.sqrt(n)
    ni = 1./max(1.e-10, n)
    return [ v[0]*ni, v[1]*ni, v[2]*ni ]

# square distance
def squareDist(p1, p2):
    """Square distance between two points."""
    dx = p1[0]-p2[0]
    dy = p1[1]-p2[1]
    dz = p1[2]-p2[2]
    return dx*dx+dy*dy+dz*dz

def dist(p1, p2):
    """Distance between two points."""
    dx = p1[0]-p2[0]
    dy = p1[1]-p2[1]
    dz = p1[2]-p2[2]
    return math.sqrt(dx*dx+dy*dy+dz*dz)

# the matrix is given by lines
def matprod(m, v):
    """Matrix-vector product."""
    [l1,l2,l3] = m # get each line of matrix
    vx = l1[0]*v[0]+l1[1]*v[1]+l1[2]*v[2]
    vy = l2[0]*v[0]+l2[1]*v[1]+l2[2]*v[2]
    vz = l3[0]*v[0]+l3[1]*v[1]+l3[2]*v[2]
    return [ vx, vy, vz ]
