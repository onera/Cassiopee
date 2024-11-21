"""Do something cool with pyTrees.
"""
import Aero
import aero
import numpy as np
try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
except:
    raise ImportError("Aero.PyTree: requires Converter.PyTree module.")
__version__ = Aero.__version__

def getTangent(z):
    """
    Makes the tangent of a 1D curve, as a new zone. The input argument shall be
    a zone. Each node of the output represents the unitary tangent vector, 
    pointing towards the tangent direction of the input 1D curve.
    """
    coords = C.getFields(Internal.__GridCoordinates__, z)
    a = Aero.getTangent(coords[0])
    return C.convertArrays2ZoneNode('line', [a])


def getDistribution(t):
    """ Given a curve, creates a 1D line of unity length with the curve's
    equivalent distribution of points. Useful when combined with Generator's
    mapping functions.
    Usage: getDistribution(t) """
    coords = C.getFields(Internal.__GridCoordinates__, t)
    a = Aero.getDistribution(coords[0])
    return C.convertArrays2ZoneNode('line', [a])


def airfoil(designation='NACA0012',Ntop=None, Nbot=None, ChordLength=1., TopDistribution=None,
            BottomDistribution=None,Closed=True,LeadingEdgePos=None):
    """Creates a 4-digit or 5-digit series NACA airfoil with discretization 
    parameters:
    designation : String of 4 or 5 digits, is the NACA airfoil identifier.    
    Ntop : Number of points of the Top side of the airfoil. 
    Nbot : Number of points of the Bottom side of the airfoil. If Nbot is not
           provided, then Ntop is the total number of points of the whole foil.
    ChordLength : The chord length of the airfoil.
    TopDistribution : A distribution dictionnary establishing the discretization
                      law of the top side of the airfoil.
    BottomDistribution : A distribution dictionnary establishing the 
                        discretization law of the bottom side of the airfoil.
    Closed : if True, it produces a closed airfoil, otherwise the curve is open.
    LeadingEdgePos : A float between 0.0 and 1.0 and typically close to 0.5
                     establishing the parametric relative position of the
                     Leading Edge position. It is used to accurately control
                     the location of the leading edge refinement point. A value
                     of 0 corresponds to the bottom side trailing edge, and a 
                     value of 1 corresponds to the top side trailing edge.
                     If None, then makes no refinement based on curve length,
                     but based only on the X coordinate.
    Usage example: airfoil(designation,Ntop,Nbot,ChordLength,TopDistribution,
                           BottomDistribution,Closed,LeadingEdgePos )
    """          
    a = Aero.airfoil(designation,Ntop,Nbot,ChordLength,TopDistribution,
                     BottomDistribution,Closed,LeadingEdgePos)
    return C.convertArrays2ZoneNode('airfoil'+designation, [a])

def linelaw(P1=(0,0,0), P2=(1,0,0), N=100, Distribution = None):
    """Create a line of N points following a given point Distribution.
    Distribution is a dictionnary with fields:
    'kind' of law: 'linear', 'tanhOneSide', 'tanhTwoSides', 'trigonometric' 
    'FirstCellHeight' -> Distance of the second point measured from P1
    'LastCellHeight' -> Distance of the before-last point measured from P2
    Usage: linelaw( (x1,y1,z1), (x2,y2,z2), N, Distribution )"""
    a = Aero.linelaw(P1, P2, N, Distribution)
    return C.convertArrays2ZoneNode('line', [a])

def splinelaw(polyLine,N=100,Distribution=None,SplineDegree=3):
    """
    Similar function to linelaw, but this time it creates a spline. polyLine
    has the points where the spline will pass through.  N is the total number of
    spline points, and Distribution is the Dictionary describing the distribution
    law, similar to linelaw, but CellHeight is not interpreted as an absolute
    value like linelaw, but rather as a non-dimensional 0 to 1 ratio of the 
    polyLine length.
    """
    polyLineA = C.getFields(Internal.__GridCoordinates__, polyLine)[0]
    a = Aero.splinelaw(polyLineA,N,Distribution,SplineDegree)
    return C.convertArrays2ZoneNode('splinelaw', [a])    

def printSmth(string):
    return Aero.printSmth(string)

# z must be a zone
def pyTreeExample(z):
    return aero.pyTreeExample(z)

# t must be a pyTree
def pyTreeExample1(t):
    return aero.pyTreeExample1(t)

def wing(sections=[airfoil(),airfoil()], span=[1.], washout=[0.], sweep=[0.], 
        dihedral=[0.], distribution=[{'points':10}],sectionShapeLaw='linear'):
    '''
    This function builds an airplane's wing section by section, using typical
    aerodynamic parameters. The result is a structured surface.

    sections : List of sections defining each airfoil, with its corresponding
    chord length. The minimum amount of sections needed to define a wing is 2
    (root and tip). Use the airfoil() function to define and discretize the 
    airfoils as you wish. Beware that the total number of points of each airfoil
    must be the same.

    span : It is a list defining the wingspan between two consecutive sections.

    washout : angle of the downwards twist of the wing at each section. You may
    only choose to give one angle of washout, and it will act as a rotation of
    the whole wing.

    sweep : backwards sweepangle of each section. Only one value imposes a single
    sweep angle for the whole wing.

    dihedral : angle of the upwards dihedral of the wing at each section. Only 
    one value imposes a single dihedral angle for the whole wing.

    distribution : Dictionnary of distribution similar to the ones used in
    linelaw. It is used to control de discretisation in the span-coordinate
    between two consecutive wing sections. If only one distribution is provided,
    the distribution will be the same for each wing element (portion of wing
    between two consecutive sections). You shall provide the keyword 'points' in
    order to impose the number of discretization points.

    sectionShapeLaw : It controls the law imposed to interpolate the airfoils 
    between two sections. Can be ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, 
    ‘quadratic‘ or ‘cubic’. Recommended value is 'linear' for trapezoidal-like
    wings. For complex geometries with non-linear interpolations, the user may
    prefer to use the function sweepSections().
    ______________
    Important note: wing supposes that the input sections are positioned in the
    XY plane. The reference point from which the washout, dihedral and sweep 
    angles are applied is the origin (0,0,0), which corresponds by default to 
    the Leading Edge of the airfoil. If you wish to re-position this reference,
    you may apply a translation to the wing sections beforehand.
    '''        
    sectionsA = []

    for i in range(len(sections)):
        sectionsA.append(C.getFields(Internal.__GridCoordinates__, sections[i])[0])
    a = Aero.wing(sectionsA, span, washout, sweep, 
        dihedral, distribution,sectionShapeLaw)
    return C.convertArrays2ZoneNode('wing', [a])   

def sweepSections(sections=[airfoil(),airfoil()], SpanPositions=None,
                  rotation=[0.], rotationLaw='linear',
                  NormalDirection=np.array([1,0,0]),
                  spine=linelaw(P2=(0,0,1)), sectionShapeLaw='linear'):     
    '''
    This function builds a sweep surface from a given profile, or a set of
    profiles, throughout a spine. The result is a structured surface.

    sections : List of sections that will be used to sweep around the spine.
    Beware that the total number of points of each section must be the same.

    SpanPositions : Vector between 0 and 1, with the same number of elements as
    the number of sections. It is used to place each section at the wished 
    position along the spine. If not provided, sections are distributed
    uniformely.

    rotation : Vector of twist angles imposed to the profiles along the spine. 
    You may use a single value, and the imposed rotation will be the same 
    throughout the spine. You may use only two values, and the rotation will be
    distributed linearly from start to end. You may use the same number of values
    as the number of sections, then the twist angle is imposed following the law
    (rotationLaw) at each section's position.

    rotationLaw : It controls the law imposed to interpolate the rotation 
    between two sections. Can be ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, 
    ‘quadratic‘ or ‘cubic’. 

    NormalDirection : 3 element vector describing the extrusion orientation
    used to find the local reference frame. As the sections are supposed to be
    placed in the XY plane, this vector will typically be (0,0,1).

    spine : 1D curve where the sections will be swept along. The spine-wise 
    discretization is exactly the same as the node distribution of the spine. 
    Hence, if you wish to control the discretization along the spine, change 
    your spine preferrably with linelaw or spinelaw.

    sectionShapeLaw : It controls the law imposed to interpolate the airfoils 
    between two sections. Can be ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, 
    ‘quadratic‘ or ‘cubic’.
    ______________
    Important note: We suppose that the input sections are positioned in the
    XY plane. The reference point from which the rotation 
    angles are applied is the origin (0,0,0). If you wish to re-position
    this reference, you may apply a translation to the sections beforehand.
    '''                        
    sectionsA = []
    for i in range(len(sections)):
        sectionsA.append(C.getFields(Internal.__GridCoordinates__, sections[i])[0])
    spineA = C.getFields(Internal.__GridCoordinates__, spine)[0]
    a = Aero.sweepSections(sectionsA, SpanPositions, rotation, rotationLaw,
                NormalDirection, spineA, sectionShapeLaw)
    return C.convertArrays2ZoneNode('sweep', [a])                        
