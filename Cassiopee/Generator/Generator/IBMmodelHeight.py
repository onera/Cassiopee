"""Compute values for F42 IBM """
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy
import math

EPSCART = 1.e-6

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## COMPUTSE INFO FOR F42 (e.g. Yplus & modelisation height etc.)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=============================================================================
# Compute the skin friction coefficient for a given emperical law
#=============================================================================
def compute_Cf(Re, Cf_law='ANSYS'):
    if Cf_law == 'ANSYS':
        return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        return 1.328*Re**(-0.5)


#=============================================================================
# Compute the corresponding yplus of a given modeling height
#=============================================================================    
def computeYplus(Re, Cf_law='ANSYS', height=0.1, L=1.):
    h0 = (L*numpy.sqrt(2))/(Re*numpy.sqrt(compute_Cf(Re,Cf_law)))
    return height/h0


#=============================================================================
# Compute the modeling height
#=============================================================================
def computeModelisationHeight(Re, Cf_law='ANSYS', yplus=100., L=1.):
    return (yplus*L*numpy.sqrt(2))/(Re*numpy.sqrt(compute_Cf(Re,Cf_law)))


#=============================================================================
# Compute the best modeling height for a given snear
#=============================================================================
def computeBestModelisationHeight(Re, h, Cf_law='ANSYS', L=1., q=1.2):
    h0 = (L*numpy.sqrt(2))/(Re*numpy.sqrt(compute_Cf(Re,Cf_law)))
    hmod = (h0-q*h)/(1.-q)
    return hmod, hmod/h0


def computeYplusOpt(Re=None,tb=None,Lref=1.,q=1.2,snear=None,Cf_law='ANSYS'):
    fail=0
    if Re is None:
        if tb is not None:
            Re = Internal.getNodeFromName(tb,"Reynolds")
            if Re is None: fail=1
            else:
                Re = Internal.getValue(Re)
        else: fail = 1
    if fail: 
        raise ValueError("computeYplusOpt: requires Reynolds number as a float or in tb.")
    fail = 0
    if snear is None:
        snear = Internal.getNodeFromName(tb,"snear")
        if snear is None: fail=1
        else: snear = Internal.getValue(snear)
    if fail:
        raise ValueError("computeYlusOpt: requires snear as a float or in tb.")

    print("Warning: estimation of the optimum y+ at Reynolds number ", Re, " and snear target at image point ", snear)
    h0 = (1.*Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law))) #Taille de maille pour y+1
    h_opti = (h0-q*snear)/(1.-q) #Hauteur de modelisation opti
    yplus_opti = h_opti/h0 #yplus opti
    # print('\nInformation for the body-fitted mesh :')
    # print('h_opti     = %1.2e'%(h_opti))
    # print('h0         = %1.2e\n'%(h0))
    # print('Information for the Cartesian mesh :')
    # print('yplus_opti = %d\n'%(int(math.ceil(yplus_opti))))
    return yplus_opti


# compute the near wall spacing in agreement with the yplus target at image points - front42
def computeSnearOpt(Re=None,tb=None,Lref=1.,q=1.2,yplus=300.,Cf_law='ANSYS'):
    fail=0
    if Re is None:
        if tb is not None:
            Re = Internal.getNodeFromName(tb,"Reynolds")
            if Re is None: fail=1
            else: Re = Internal.getValue(Re)
        else: fail = 1
    if fail: 
        raise ValueError("computeSnearOpt: requires Reynolds number as a float or in tb.")


    print("Estimation of the optimum near-wall spacing at Reynolds number ", Re, " and yplus target at image point ", yplus)
    h_mod = (yplus*Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law)))
    h0    = (Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law))) #Taille de maille pour y+=1
    n     = int(math.ceil(math.log(1-yplus*(1-q))/math.log(q))) # number of cells in the BF mesh for the height h
    snear_opti = q**(n-1)*h0 # best snear for the target yplus
    # print('\nInformation for the body-fitted mesh :')
    # print('h           = %1.2e'%(h_mod))
    # print('h0          = %1.2e\n'%(h0))
    # print('Information for the Cartesian mesh :')
    # print('snear_opti  = %1.3e\n'%(snear_opti))
    return snear_opti


def getMinimumCartesianSpacing(t):
    baseC = Internal.getNodeFromName1(t, 'CARTESIAN')
    if baseC is None: return -1.

    zonesC = Internal.getZones(baseC)
    dxmin = 1.e6
    for z in zonesC:
        dx = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
        if dx < dxmin: dxmin = dx

    print('Minimum spacing on Cartesian grids = %f.'%dxmin)
    return dxmin
