"""Compute values for F42 IBM """
import Converter.PyTree as C
import Converter.Internal as Internal
import Geom.IBM as D_IBM

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## COMPUTE INFO FOR F42 (e.g. Yplus & modelisation height etc.)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=============================================================================
# Compute the skin friction coefficient for a given emperical law
#=============================================================================
def compute_Cf(Re, Cf_law='ANSYS'):
    return D_IBM.compute_Cf(Re, Cf_law=Cf_law)

#=============================================================================
# Compute the corresponding yplus of a given modeling height
#=============================================================================
def computeYplus(Re, Cf_law='ANSYS', height=0.1, L=1.):
    return D_IBM.computeYplus(Re, Cf_law=Cf_law, height=height, L=L)

#=============================================================================
# Compute the modeling height
#=============================================================================
def computeModelisationHeight(Re, Cf_law='ANSYS', yplus=100., L=1.):
    return D_IBM.computeModelisationHeight(Re, Cf_law=Cf_law, yplus=yplus, L=L)

#=============================================================================
# Compute the best modeling height for a given snear
#=============================================================================
def computeBestModelisationHeight(Re, h, Cf_law='ANSYS', L=1., q=1.2):
    return D_IBM.computeBestModelisationHeight(Re, h, Cf_law=Cf_law, L=L, q=q)

def computeYplusOpt(Re=None,tb=None,Lref=1.,q=1.2,snear=None,Cf_law='ANSYS'):
    return D_IBM.computeYplusOpt(Re=Re,tb=tb,Lref=Lref,q=q,snear=snear,Cf_law=Cf_law)

# compute the near wall spacing in agreement with the yplus target at image points - front42
def computeSnearOpt(Re=None,tb=None,Lref=1.,q=1.2,yplus=300.,Cf_law='ANSYS'):
    return computeSnearOpt(Re=Re,tb=tb,Lref=Lref,q=q,yplus=yplus,Cf_law=Cf_law)

def getMinimumCartesianSpacing(t):
    return D_IBM.getMinimumCartesianSpacing(t)