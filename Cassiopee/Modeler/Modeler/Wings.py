"""Do something cool with arrays.
"""
__version__ = '2.2'
__author__ = "You"

import numpy
import Converter.Array3D

try: import Generator as G
except: raise ImportError("Aero: requires Generator module.")
try: import Geom as D
except: raise ImportError("Aero: requires Geom module.")
try: import scipy.interpolate
except: raise ImportError("Aero: requires scipy interpolate module.")
try: import Transform as T
except: raise ImportError("Aero: requires Transform module.")


def airfoil(designation='NACA0012', Ntop=None, Nbot=None, ChordLength=1.,
            TopDistribution=None, BottomDistribution=None,
            Closed=True, LeadingEdgePos=None):
    """Creates a 4-digit or 5-digit series NACA airfoil with discretization
    parameters:
    designation: String of 4 or 5 digits, is the NACA airfoil identifier.
    Ntop: Number of points of the Top side of the airfoil.
    Nbot: Number of points of the Bottom side of the airfoil. If Nbot is not
          provided, then Ntop is the total number of points of the whole foil.
    ChordLength: The chord length of the airfoil.
    TopDistribution: A distribution dictionnary establishing the discretization
                     law of the top side of the airfoil.
    BottomDistribution: A distribution dictionnary establishing the
                        discretization law of the bottom side of the airfoil.
    Closed: if True, it produces a closed airfoil, otherwise the curve is open.
    LeadingEdgePos: A float between 0.0 and 1.0 and typically close to 0.5
                    establishing the parametric relative position of the
                    Leading Edge position. It is used to accurately control
                  the location of the leading edge refinement point. A value
                  of 0 corresponds to the bottom side trailing edge, and a
                  value of 1 corresponds to the top side trailing edge.
                  If None, then makes no refinement based on curve length,
                  but based only on the X coordinate.
    Usage example: airfoil(designation, Ntop, Nbot, ChordLength, TopDistribution,
                           BottomDistribution, Closed, LeadingEdgePos )
    """

    # Prepares the X-distributions
    if not Ntop and not Nbot and not designation.find('.') != -1:
        Ntop = 50
        Nbot = 50
    elif not not Ntop and not Nbot:
        Nbot = Ntop/2 + Ntop%2
        Ntop /= 2

    if not TopDistribution:
        TopDistribution = {'kind':'trigonometric','parameter':2}
    if not BottomDistribution:
        if TopDistribution['kind'] == 'tanhTwoSides':
            BottomDistribution = TopDistribution.copy()
            BottomDistribution['FirstCellHeight'] = TopDistribution['LastCellHeight']
            BottomDistribution['LastCellHeight'] = TopDistribution['FirstCellHeight']
        elif TopDistribution['kind'] == 'tanhOneSide':
            BottomDistribution = TopDistribution.copy()
        elif TopDistribution['kind'] == 'trigonometric':
            BottomDistribution = TopDistribution.copy()
            BottomDistribution['parameter'] *= -1

    NACAstringLoc = designation.find('NACA')
    # Determines the kind of airfoil to generate
    if designation.find('.') != -1: # Then user wants to import an airfoil from file
        Imported = numpy.genfromtxt(designation, dtype=numpy.float, skip_header=0, usecols=(0,1))
        # Deletes useless lines
        RowsToDelete = []
        for i in range(len(Imported[:,0])):
            if any(numpy.isnan(Imported[i])) or any(Imported[i]>1.5):
                RowsToDelete.append(i)
        Imported = numpy.delete(Imported, RowsToDelete, axis=0)
        # Checks for the format of the coordinate points
        Monot = numpy.diff(Imported[:,0])
        MonotIND = numpy.where(Monot < 0)[0]
        if len(MonotIND) == 1:
            # Lednicer format: Both sides start from Leading Edge
            if Imported[MonotIND,1] < Imported[-1,1]:
                xL = numpy.flipud(Imported[:MonotIND+1,0])
                yL = numpy.flipud(Imported[:MonotIND+1,1])
                xU = Imported[MonotIND+1:,0]
                yU = Imported[MonotIND+1:,1]

            else:
                xU = Imported[:MonotIND+1,0]
                yU = Imported[:MonotIND+1,1]
                xL = numpy.flipud(Imported[MonotIND+1:,0])
                yL = numpy.flipud(Imported[MonotIND+1:,1])
        else:
            # Selig format: Starts and ends from trailing edge
            if Imported[1,1] > Imported[-2,1]: Imported = numpy.flipud(Imported)
            xMin = numpy.argmin(Imported[:,0])
            xL= Imported[:xMin+1,0]
            yL = Imported[:xMin+1,1]
            xU= Imported[xMin:,0]
            yU = Imported[xMin:,1]

        if not not LeadingEdgePos and not not Ntop: # then it makes a dense spline sampling.
            # Therefore, for fine refinement the user must provide both Ntop and
            # LeadingEdgePos. Refinement is only made through curvilinear abcissa.
            tU = numpy.zeros(xU.shape)
            tL = numpy.zeros(xL.shape)
            tU[1:] = numpy.sqrt((numpy.diff(xU))**2 + (numpy.diff(yU))**2)
            tL[1:] = numpy.sqrt((numpy.diff(xL))**2 + (numpy.diff(yL))**2)
            tU = numpy.cumsum(tU)
            tL = numpy.cumsum(tL)
            tU /= tU[-1]
            tL /= tL[-1]
            nt = numpy.linspace(0, 1, 5000)
            Degree = 3
            xL = scipy.interpolate.spline(tL, xL, nt, Degree)
            yL = scipy.interpolate.spline(tL, yL, nt, Degree)
            xU = scipy.interpolate.spline(tU, xU, nt, Degree)
            yU = scipy.interpolate.spline(tU, yU, nt, Degree)

        if Closed: yL[0]=yU[-1]=0.5*(yL[0]+yU[-1])
        Airfoil = D.line((0,0,0), (1,0,0), len(xL)+len(xU)-1 )
        Airfoil[1][0] = numpy.hstack((xL,xU[1:]))
        Airfoil[1][1] = numpy.hstack((yL,yU[1:]))
    elif NACAstringLoc != -1: # Then user wants to generate a naca-series airfoil
        if not Ntop: Ntop = Nbot = 100
        xU = linelaw((0,0,0), (1,0,0), Ntop, TopDistribution)[1][0]
        xL = linelaw((1,0,0), (0,0,0), Nbot, BottomDistribution)[1][0]
        NACAidentifier = designation[NACAstringLoc + 4:len(designation)]
        # NACA constants
        a0= 0.2969; a1=-0.1260; a2=-0.3516; a3= 0.2843; a4=-0.1036 if Closed else -0.1015

        if len(NACAidentifier) == 4: # 4-digit NACA
            m = float(NACAidentifier[0])*0.01  # Maximum camber
            p = float(NACAidentifier[1])*0.1   # Maximum camber location
            t = float(NACAidentifier[2:])*0.01 # Maximum thickness
            ytU = 5.*t*(a0*numpy.sqrt(xU)+a1*(xU)+a2*(xU)**2+a3*(xU)**3+a4*(xU)**4)
            ytL = 5.*t*(a0*numpy.sqrt(xL)+a1*(xL)+a2*(xL)**2+a3*(xL)**3+a4*(xL)**4)
            if m == 0: # no cambered airfoil, it is symmetric
                if Closed: ytU[-1]=ytL[0]=0 # Forces the closure of the airfoil at trailing edge
                Airfoil = D.line((0,0,0), (1,0,0), Ntop+Nbot-1)
                Airfoil[1][0] = numpy.hstack((xL,xU[1:]))
                Airfoil[1][1] = numpy.hstack((-ytL,ytU[1:]))
            else:      # cambered airfoil, non-symmetric
                ycU = numpy.zeros(Ntop)
                ycL = numpy.zeros(Nbot)
                ycU[xU<=p]= m*( xU[xU<=p]/(p**2) )*(2.*p-(xU[xU<=p]))
                ycU[xU>p]= m*( (1.-xU[xU>p])/((1.-p)**2) )*(1.-2.*p+(xU[xU>p]))
                ycL[xL<=p]= m*( xL[xL<=p]/(p**2) )*(2.*p-(xL[xL<=p]))
                ycL[xL>p]= m*( (1.-xL[xL>p])/((1.-p)**2) )*(1.-2.*p+(xL[xL>p]))
                thU = numpy.zeros(Ntop)
                thL = numpy.zeros(Nbot)
                thU[xU<=p]= (2.*m/(p**2))*(p-(xU[xU<=p]))
                thU[xU>p]= (2.*m/((1.-p)**2))*(p-(xU[xU>p]))
                thL[xL<=p]= (2.*m/(p**2))*(p-(xL[xL<=p]))
                thL[xL>p]= (2.*m/((1.-p)**2))*(p-(xL[xL>p]))
                thU = numpy.arctan(thU); thL = numpy.arctan(thL)
                xU = xU - ytU*numpy.sin(thU); yU = ycU + ytU*numpy.cos(thU)
                xL = xL + ytL*numpy.sin(thL); yL = ycL - ytL*numpy.cos(thL)
                if Closed: yU[-1]=yL[0] # Forces the closure of the airfoil at trailing edge
                Airfoil = D.line((0,0,0), (1,0,0), Ntop+Nbot-1)
                Airfoil[1][0] = numpy.hstack((xL,xU[1:]))
                Airfoil[1][1] = numpy.hstack((yL,yU[1:]))

        elif len(NACAidentifier) == 5: # 5-digit NACA
            cld = float(NACAidentifier[0]) *(3./2.)*0.1
            p = float(NACAidentifier[1])
            if p > 5:
                print ('Warning: second digit of 5-digit NACA identifier > 5, switched to 5')
                p = 5
            p   /= 20.
            q   = int(NACAidentifier[2])
            t   = float(NACAidentifier[3:])*0.01
            if q == 0: # standard
                P   = numpy.array([  0.05,     0.1,     0.15,    0.2,     0.25  ])
                R   = numpy.array([  0.0580,   0.1260,  0.2025,  0.2900,  0.3910])
                K   = numpy.array([361.4,     51.64,   15.957,   6.643,   3.230 ])
            else: # reflex
                P   = numpy.array([  0.1,      0.15,    0.2,     0.25  ])
                R   = numpy.array([  0.13,     0.217,   0.318,   0.441 ])
                K   = numpy.array([ 51.99,    15.793,   6.520,   3.191 ])
                K2K1= numpy.array([  0.000764, 0.00677, 0.0303,  0.1355])

            ytU = 5.*t*(a0*numpy.sqrt(xU)+a1*(xU)+a2*(xU)**2+a3*(xU)**3+a4*(xU)**4)
            ytL = 5.*t*(a0*numpy.sqrt(xL)+a1*(xL)+a2*(xL)**2+a3*(xL)**3+a4*(xL)**4)
            if p == 0: # no cambered airfoil, it is symmetric
                if Closed: ytU[-1]=ytL[0] # Forces the closure of the airfoil at trailing edge
                Airfoil = D.line((0,0,0), (1,0,0), Ntop+Nbot-1)
                Airfoil[1][0] = numpy.hstack((xL,xU[1:]))
                Airfoil[1][1] = numpy.hstack((-ytL,ytU[1:]))
            else:      # cambered airfoil, non-symmetric

                inter_pr = scipy.interpolate.UnivariateSpline(P,R)
                inter_pk = scipy.interpolate.UnivariateSpline(P,K)
                r = inter_pr(p)
                k1= inter_pk(p)
                ycU = numpy.zeros(Ntop)
                ycL = numpy.zeros(Nbot)
                Scale = cld/0.3
                thU = numpy.zeros(Ntop)
                thL = numpy.zeros(Nbot)
                if q==0: # standard equations
                    ycU[xU<=r]= Scale*(k1/6.)*((xU[xU<=r])**3 - 3*r*(xU[xU<=r])**2 + (r**2)*(3-r)*(xU[xU<=r]) )
                    ycU[xU>r]= Scale*(k1/6.)*(r**3)*(1-(xU[xU>r]))
                    ycL[xL<=r]= Scale*(k1/6.)*((xL[xL<=r])**3 - 3*r*(xL[xL<=r])**2 + (r**2)*(3-r)*(xL[xL<=r]) )
                    ycL[xL>r]= Scale*(k1/6.)*(r**3)*(1-(xL[xL>r]))
                    thU[xU<=r]= Scale*(k1/6.)* ( 3.*(xU[xU<=r])**2 - 6.*r*(xU[xU<=r]) + (r**2)*(3.-r) )
                    thU[xU>r]= -Scale*(k1/6.)*(r**3)*(xU[xU>r]*0)
                    thL[xL<=r]= Scale*(k1/6.)* ( 3.*(xL[xL<=r])**2 - 6.*r*(xL[xL<=r]) + (r**2)*(3.-r) )
                    thL[xL>r]= -Scale*(k1/6.)*(r**3)*(xL[xL>r]*0)
                else:   # reflex equations
                    inter_pk2k1 = scipy.interpolate.UnivariateSpline(P,K2K1)
                    k2k1 = inter_pk2k1(p)
                    ycU[xU<=r]= Scale*(k1/6.)*((xU[xU<=r] - r)**3 -k2k1*(xU[xU<=r])*((1-r)**3 - r**3 ) + r**3 )
                    ycU[xU>r]= Scale*(k1/6.)*(k2k1*(xU[xU>r] - r)**3 -k2k1*(xU[xU>r])*((1-r)**3 - r**3 ) + r**3 )
                    ycL[xL<=r]= Scale*(k1/6.)*((xL[xL<=r] - r)**3 -k2k1*(xL[xL<=r])*((1-r)**3 - r**3 ) + r**3 )
                    ycL[xL>r]= Scale*(k1/6.)*(k2k1*(xL[xL>r] - r)**3 -k2k1*(xL[xL>r])*((1-r)**3 - r**3 ) + r**3 )
                    thU[xU<=r]= Scale*(k1/6.)*(3.*(xU[xU<=r] - r)**2 -k2k1*((1-r)**3) - r**3 )
                    thU[xU>r]= Scale*(k1/6.)*(k2k1*3.*(xU[xU>r] - r)**2 -k2k1*((1-r)**3) - r**3 )
                    thL[xL<=r]= Scale*(k1/6.)*(3.*(xL[xL<=r] - r)**2 -k2k1*((1-r)**3) - r**3 )
                    thL[xL>r]= Scale*(k1/6.)*(k2k1*3.*(xL[xL>r] - r)**2 -k2k1*((1-r)**3) - r**3 )
                thU = numpy.arctan(thU); thL = numpy.arctan(thL)
                xU = xU - ytU*numpy.sin(thU); yU = ycU + ytU*numpy.cos(thU)
                xL = xL + ytL*numpy.sin(thL); yL = ycL - ytL*numpy.cos(thL)
                if Closed: yU[-1]=yL[0] # Forces the closure of the airfoil at trailing edge
                Airfoil = D.line((0,0,0), (1,0,0), Ntop+Nbot-1)
                Airfoil[1][0] = numpy.hstack((xL,xU[1:]))
                Airfoil[1][1] = numpy.hstack((yL,yU[1:]))
    else:
        print ('airfoil: designation "%s" not recognized.'%designation)
        return -1
    # Scaling
    Airfoil[1] *= ChordLength
    if not LeadingEdgePos or not Ntop: # Airfoil with standard X-based refinement
        return Airfoil
    else: # Final refinement around airfoil's curve length
        DistrCurve = D.line((0,0,0), (1,0,0), Ntop+Nbot-1)
        if (type(TopDistribution) is dict) and (type(BottomDistribution) is dict):
            keysTD = TopDistribution.keys()
            keysBD = BottomDistribution.keys()
            if 'FirstCellHeight' in keysTD: TopDistribution['FirstCellHeight'] *= 1-LeadingEdgePos
            if 'LastCellHeight' in keysTD: TopDistribution['LastCellHeight'] *= 1-LeadingEdgePos
            if 'FirstCellHeight' in keysBD: BottomDistribution['FirstCellHeight'] *= LeadingEdgePos
            if 'LastCellHeight' in keysBD: BottomDistribution['LastCellHeight'] *= LeadingEdgePos
            TopCurveDist = linelaw((LeadingEdgePos,0,0),(1,0,0),Ntop,TopDistribution)
            BotCurveDist = linelaw((0,0,0),(LeadingEdgePos,0,0),Nbot,BottomDistribution)
            DistrCurve[1][0] = numpy.hstack((BotCurveDist[1][0],TopCurveDist[1][0][1:]))
        FoilDistribution = D.getDistribution(DistrCurve)
        return G.map1d(Airfoil,FoilDistribution)


def linelaw(P1=(0,0,0), P2=(1,0,0), N=100, Distribution=None):
    """Create a line of N points following a given point Distribution.
    Distribution is a dictionnary with fields:
    'kind' of law: 'linear', 'tanhOneSide', 'tanhTwoSides', 'trigonometric'
    'FirstCellHeight' -> Distance of the second point measured from P1
    'LastCellHeight' -> Distance of the before-last point measured from P2
    Usage: linelaw( (x1,y1,z1), (x2,y2,z2), N, Distribution )"""
    if not Distribution: return D.line(P1, P2, N)
    else:
        Line = D.line(P1,P2,N)
        if Distribution['kind'] == 'tanhOneSide':
            Length = numpy.linalg.norm(numpy.array(P2)-numpy.array(P1))
            dy = Distribution['FirstCellHeight']/Length
            Dir = numpy.array([P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]])/Length
            S = G.getTanhDist__(N,dy,dy,1)[1][0]*Length
            Height = S[1]-S[0]
            ErrorHeight = abs(100*(1-Height/(dy*Length)))
            if ErrorHeight > 1.:
                print ('--------')
                print ('Warning: Distribution of kind tanhOneSide resulted in an')
                print ('effective cell Height of:',Height,' which differs from the')
                print ('desired one,',dy*Length,' a relative amount of:',ErrorHeight,'%.')
                print ('Try different discretization parameters for better result.')
                print ('--------')

            Line[1][0] = S*Dir[0]+P1[0]
            Line[1][1] = S*Dir[1]+P1[1]
            Line[1][2] = S*Dir[2]+P1[2]
            # Strictly represents boundaries (no epsilon deviation):
            Line[1][0][0] = P1[0]
            Line[1][1][0] = P1[1]
            Line[1][2][0] = P1[2]
            Line[1][0][-1] = P2[0]
            Line[1][1][-1] = P2[1]
            Line[1][2][-1] = P2[2]
            return Line
        elif Distribution['kind'] == 'tanhTwoSides':
            Length = numpy.linalg.norm(numpy.array(P2)-numpy.array(P1))
            dy = [0.,0.]
            dy[0] = Distribution['FirstCellHeight']/Length
            dy[1] = Distribution['LastCellHeight']/Length
            Dir = numpy.array([P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]])/Length
            S = G.getTanhDist__(N,dy[0],dy[1],2)[1][0]*Length
            Height1 = S[1]-S[0]; Height2 = S[-1]-S[-2]
            ErrorHeight1 = abs(100*(1-Height1/(dy[0]*Length)))
            ErrorHeight2 = abs(100*(1-Height2/(dy[1]*Length)))
            if ErrorHeight1 > 1.:
                print ('--------')
                print ('Warning: Distribution of kind tanhTwoSides resulted in an')
                print ('effective first cell Height of:',Height1,' which differs from the')
                print ('desired one,',dy[0]*Length,' a relative amount of:',ErrorHeight1,'%.')
                print ('Try different discretization parameters for better result.')
                print ('--------')
            elif ErrorHeight2 > 1.:
                print ('--------')
                print ('Warning: Distribution of kind tanhTwoSides resulted in an')
                print ('effective last cell Height of:',Height2,' which differs from the')
                print ('desired one,',dy[1]*Length,' a relative amount of:',ErrorHeight2,'%.')
                print ('Try different discretization parameters for better result.')
                print ('--------')
            Line[1][0] = S*Dir[0]+P1[0]
            Line[1][1] = S*Dir[1]+P1[1]
            Line[1][2] = S*Dir[2]+P1[2]
            # Strictly represents boundaries (no epsilon deviation):
            Line[1][0][0] = P1[0]
            Line[1][1][0] = P1[1]
            Line[1][2][0] = P1[2]
            Line[1][0][-1] = P2[0]
            Line[1][1][-1] = P2[1]
            Line[1][2][-1] = P2[2]
            return Line
        elif Distribution['kind'] == 'trigonometric':
            Length = numpy.linalg.norm(numpy.array(P2)-numpy.array(P1))
            p = Distribution['parameter']
            Dir = numpy.array([P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]])/Length
            S = getTrigoLinDistribution__(N, p)*Length
            Line[1][0] = S*Dir[0]+P1[0]
            Line[1][1] = S*Dir[1]+P1[1]
            Line[1][2] = S*Dir[2]+P1[2]
            # Strictly represents boundaries (no epsilon deviation):
            Line[1][0][0] = P1[0]
            Line[1][1][0] = P1[1]
            Line[1][2][0] = P1[2]
            Line[1][0][-1] = P2[0]
            Line[1][1][-1] = P2[1]
            Line[1][2][-1] = P2[2]
            return Line
        else:
            return Line


def getTrigoLinDistribution__(Nx, p):
    """Private function used in line. Returns a [0,1] vector of Nx points
    following a mix of trigonometric and linear distribution. The parameter
    p in [-3,3] controls the discretization law."""

    x = numpy.linspace(0,1,Nx)
    if p <= 3 and p >= 2:
        L = x
        S1 = 1. + numpy.sin((numpy.pi/2.)*(x-1.))
        a = 3 - p
        b = p - 2
        return b*L + a*S1
    elif p <= 2 and p >= 1:
        S1 = 1. + numpy.sin((numpy.pi/2.)*(x-1.))
        C = 0.5*(1.+numpy.cos(numpy.pi*(x+1)))
        a = 2 - p
        b = p - 1
        return b*S1 + a*C
    elif p <= 1 and p >= 0:
        C = 0.5*(1.+numpy.cos(numpy.pi*(x+1)))
        L = x
        a = 1 - p
        b = p
        return b*C + a*L
    elif p <= 0 and p >= -1:
        L = x
        C = 0.5*(1.+numpy.cos(numpy.pi*(x+1)))
        a = - p
        b = p - (-1)
        return b*L + a*C
    elif p <= -1 and p >= -2:
        C = 0.5*(1.+numpy.cos(numpy.pi*(x+1)))
        S2 = numpy.sin((numpy.pi/2.)*x)
        a = (-1) - p
        b = p - (-2)
        return b*C + a*S2
    elif p <= -2 and p >= -3:
        S2 = numpy.sin((numpy.pi/2.)*x)
        L = x
        a = (-2) - p
        b = p - (-3)
        return b*S2 + a*L
    else:
        print ('Warning: parameter p=',p,'out of allowed bounds [3,-3]. Switched to p=0.')
        return x


def splinelaw(polyLine, N=100, Distribution=None, SplineDegree=3):
    """
    Similar function to linelaw, but this time it creates a spline. polyLine
    has the points where the spline will pass through.  N is the total number of
    spline points, and Distribution is the Dictionary describing the distribution
    law, similar to linelaw, but CellHeight is not interpreted as an absolute
    value like linelaw, but rather as a non-dimensional 0 to 1 ratio of the
    polyLine length.
    """
    shPL = polyLine[1].shape
    if shPL[1] >= N:
        raise ValueError("Aero splinelaw: spline point number (%d) is too low \
        for the polyline provided (%d)"%(N,shPL[1]))
    s = numpy.zeros(shPL)
    s[:,1:] += numpy.linalg.norm(numpy.diff(polyLine[1],axis=1),axis=0)
    s = numpy.cumsum(s,axis=1)
    s /= s[:,-1][numpy.newaxis].T
    if not Distribution:
        nt = numpy.linspace(0., 1., N)
    elif 'linear' in Distribution['kind']:
        nt = numpy.linspace(0., 1., N)
    elif 'trigonometric' in Distribution['kind']:
        p = 2 if not 'parameter' in Distribution else Distribution['parameter']
        nt = getTrigoLinDistribution__(N, p)
    elif 'tanhOneSide' in Distribution['kind']:
        if not 'FirstCellHeight' in Distribution:
            print ('splinelaw: parameter FirstCellHeight not found, switching to linear.')
            nt = numpy.linspace(0., 1., N)
        else:
            CellH = Distribution['FirstCellHeight']
            if CellH > 1 or CellH < 0:
                print ('splinelaw: selected ratio FirstCellHeight of %g is out of the bounds [0,1], switch to linear.'%(CellH))
                nt = numpy.linspace(0., 1., N)
            else:
                nt = G.getTanhDist__(N,CellH,CellH,1)[1][0]
    elif 'tanhTwoSides' in Distribution['kind']:
        if not 'FirstCellHeight' in Distribution:
            print ('splinelaw: parameter FirstCellHeight not found, switching to linear.')
            nt = numpy.linspace(0., 1., N)
        else:
            CellHstart = Distribution['FirstCellHeight']
            CellHend = CellHstart if not 'LastCellHeight' in Distribution else Distribution['LastCellHeight']
            if CellHstart > 1 or CellHstart < 0:
                print ('splinelaw: selected ratio FirstCellHeight of %g is out of the bounds [0,1], switch to linear.'%(CellHstart))
                nt = numpy.linspace(0., 1., N)
            elif CellHend > 1 or CellHend < 0:
                print ('splinelaw: selected ratio LastCellHeight of %g is out of the bounds [0,1], switch to linear.'%(CellHend))
                nt = numpy.linspace(0., 1., N)
            else:
                nt = G.getTanhDist__(N,CellHstart,CellHend,2)[1][0]

    splX = scipy.interpolate.InterpolatedUnivariateSpline(s[0,:], polyLine[1][0,:], k=SplineDegree)
    splY = scipy.interpolate.InterpolatedUnivariateSpline(s[1,:], polyLine[1][1,:], k=SplineDegree)
    splZ = scipy.interpolate.InterpolatedUnivariateSpline(s[2,:], polyLine[1][2,:], k=SplineDegree)
    SplineX = splX(nt)
    SplineY = splY(nt)
    SplineZ = splZ(nt)

    return ['x,y,z', numpy.vstack((SplineX,SplineY,SplineZ)), N, 1, 1]


def wing(sections=[airfoil(),airfoil()], span=[1.], washout=[0.], sweep=[0.],
        dihedral=[0.], distribution=[{'points':10}],sectionShapeLaw='linear'):
    '''
    This function builds an airplane's wing section by section, using typical
    aerodynamic parameters. The result is a structured surface.

    sections: List of sections defining each airfoil, with its corresponding
    chord length. The minimum amount of sections needed to define a wing is 2
    (root and tip). Use the airfoil() function to define and discretize the
    airfoils as you wish. Beware that the total number of points of each airfoil
    must be the same.

    span: It is a list defining the wingspan between two consecutive sections.

    washout: angle of the downwards twist of the wing at each section. You may
    only choose to give one angle of washout, and it will act as a rotation of
    the whole wing.

    sweep: backwards sweepangle of each section. Only one value imposes a single
    sweep angle for the whole wing.

    dihedral: angle of the upwards dihedral of the wing at each section. Only
    one value imposes a single dihedral angle for the whole wing.

    distribution: Dictionnary of distribution similar to the ones used in
    linelaw. It is used to control de discretisation in the span-coordinate
    between two consecutive wing sections. If only one distribution is provided,
    the distribution will be the same for each wing element (portion of wing
    between two consecutive sections). You shall provide the keyword 'points' in
    order to impose the number of discretization points.

    sectionShapeLaw: It controls the law imposed to interpolate the airfoils
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

    # Checks the number of sections and gets the number of wing elements
    Ns = len(sections)

    WingElts = Ns-1
    if WingElts == 0:
        sections.append(sections[0])
        WingElts = 1
        Ns += 1
    if len(span) != WingElts: span = span[0]*numpy.ones(WingElts)

    # Makes the vector of sections spanwise positions
    SectionsSpanPositions = numpy.hstack((0.,numpy.cumsum(span)))
    print ('SectionsSpanPositions:',SectionsSpanPositions,'Ns=',Ns)
    # Applies the washout, sweep and dihedral angles to each section
    Dx = 0.; Dy = 0.
    for s in range(1,Ns):
        wI = 0 if len(washout)!=(Ns-1) else s-1
        dI = 0 if len(dihedral)!=(Ns-1) else s-1
        sI = 0 if len(sweep)!=(Ns-1) else s-1
        spI = 0 if len(span)!=(Ns-1) else s-1
        sectionChord = numpy.max(sections[s][1][0]) - numpy.min(sections[s][1][0])
        sections[s] = T.rotate(sections[s],(0.5*sectionChord,0,0),(0,0,washout[wI]))
        Dx += span[spI]*numpy.sqrt(1./numpy.cos(numpy.deg2rad(sweep[sI])) -1.)
        Dy += span[spI]*numpy.sqrt(1./numpy.cos(numpy.deg2rad(dihedral[dI])) -1.)
        sections[s] = T.translate(sections[s],(Dx,Dy,0))


    # Invokes the wing
    Ne = numpy.zeros(WingElts,numpy.int32)
    for w in range(WingElts):
        disI = 0 if len(distribution) != WingElts else w
        Ne[w] = distribution[disI]['points']
    Ntot = numpy.sum(Ne) - (Ns - 2)
    Wing = Converter.Array3D.convertArrays2Arrays3D(T.addkplane([sections[0]],Ntot-1))[0]

    # Makes the interpolation matrices based upon the provided sections
    W4dims = numpy.shape(Wing[1])
    print ('W4dims:',W4dims,'Ntot:',Ntot)
    InterpXmatrix = numpy.zeros((Ns,W4dims[1]),dtype=numpy.float64,order='F')
    InterpYmatrix = numpy.zeros((Ns,W4dims[1]),dtype=numpy.float64,order='F')
    for s in range(Ns):
        InterpXmatrix[s,:] = sections[s][1][0]
        InterpYmatrix[s,:] = sections[s][1][1]

    # Searches the spanwise coordinate
    ElementSpanPositions = numpy.zeros(Ntot,numpy.float64)
    k0 = 0
    PreviousSpan = 0. # and this is its corresponding spanwise position
    for w in range(WingElts):
#        print 'k0:',k0,'PreviousSpan',PreviousSpan,'k0+Ne[w]:',k0+Ne[w]
        disI = 0 if len(distribution) != WingElts else w
        Interpolation = 'linear' if not 'kind' in distribution[disI] else distribution[disI]['kind']
        if Interpolation == 'linear':
            ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*numpy.linspace(0.,1.,Ne[w])
        elif Interpolation == 'trigonometric':
            p = 2 if not 'parameter' in distribution[disI] else distribution[disI]['parameter']
            ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*getTrigoLinDistribution__(Ne[w], p)
        elif Interpolation == 'tanhOneSide':
            if not 'FirstCellHeight' in distribution[disI]:
                print ('wing: parameter FirstCellHeight not found in wing element %d, switching to linear.'%(w))
                ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*numpy.linspace(0.,1.,Ne[w])
            else:
                CellH = distribution[disI]['FirstCellHeight']/span[w]
                ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*G.getTanhDist__(Ne[w],CellH,CellH,1)[1][0]
        elif Interpolation == 'tanhTwoSides':
            if not 'FirstCellHeight' in distribution[disI]:
                print ('wing: parameter FirstCellHeight not found in wing element %d, switching to linear.'%(w))
                ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*numpy.linspace(0.,1.,Ne[w])
            else:
                CellHstart = distribution[disI]['FirstCellHeight']/span[w]
                CellHend = CellHstart if not 'LastCellHeight' in distribution[disI] else distribution[disI]['LastCellHeight']/span[w]
                ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*G.getTanhDist__(Ne[w],CellHstart,CellHend,2)[1][0]
        else:
            print ("Warning: wing distribution kind '%s' not recognized for wing element %d. Making 'linear' instead."%(Interpolation,w))
            ElementSpanPositions[k0:k0+Ne[w]] = PreviousSpan + span[w]*numpy.linspace(0.,1.,Ne[w])
        PreviousSpan = ElementSpanPositions[k0+Ne[w]-1]
        k0 += Ne[w] - 1
#    print 'ElementSpanPositions:',ElementSpanPositions

    # Makes the interpolation functions
#    print 'InterpXmatrix shape',InterpXmatrix.shape

    # Old version, makes
    interpX = scipy.interpolate.interp1d( SectionsSpanPositions, InterpXmatrix, axis=0,
        kind=sectionShapeLaw, bounds_error=False, fill_value=InterpXmatrix[0,:])
    interpY = scipy.interpolate.interp1d( SectionsSpanPositions, InterpYmatrix, axis=0,
        kind=sectionShapeLaw, bounds_error=False, fill_value=InterpYmatrix[0,:])

    # Constructs the wing section by section
    for k in range(Ntot):
        # X and Y coords are based upon the interpolating functions
#        print 'shape interpX(ElementSpanPositions[k]):',interpX(ElementSpanPositions[k]).shape
#        print 'shape Wing[1][0][:,:,k]',Wing[1][0][:,:,k].shape
        Wing[1][0][:,:,k] = interpX(ElementSpanPositions[k])[numpy.newaxis].T#numpy.reshape(interpX(ElementSpanPositions[k]).flatten(),(W4dims[1],W4dims[2]),order='F')
        Wing[1][1][:,:,k] = interpY(ElementSpanPositions[k])[numpy.newaxis].T#numpy.reshape(interpY(ElementSpanPositions[k]).flatten(),(W4dims[1],W4dims[2]),order='F')
        # Z coords:
        Wing[1][2][:,:,k] = numpy.ones((W4dims[1],W4dims[2]),order='F')*ElementSpanPositions[k]

    return Converter.Array3D.convertArrays3D2Arrays([Wing])[0]


def sweepSections(sections=[airfoil(),airfoil()], SpanPositions=None,
                  rotation=[0.], rotationLaw='linear',
                  NormalDirection=numpy.array([1,0,0]),
                  spine=linelaw(P2=(0,0,1)), sectionShapeLaw='linear'):
    '''
    This function builds a sweep surface from a given profile, or a set of
    profiles, throughout a spine. The result is a structured surface.

    sections: List of sections that will be used to sweep around the spine.
    Beware that the total number of points of each section must be the same.

    SpanPositions: Vector between 0 and 1, with the same number of elements as
    the number of sections. It is used to place each section at the wished
    position along the spine. If not provided, sections are distributed
    uniformely.

    rotation: Vector of twist angles imposed to the profiles along the spine.
    You may use a single value, and the imposed rotation will be the same
    throughout the spine. You may use only two values, and the rotation will be
    distributed linearly from start to end. You may use the same number of values
    as the number of sections, then the twist angle is imposed following the law
    (rotationLaw) at each section's position.

    rotationLaw: It controls the law imposed to interpolate the rotation
    between two sections. Can be ‘linear’, ‘nearest’, ‘zero’, ‘slinear’,
    ‘quadratic‘ or ‘cubic’.

    NormalDirection: 3 element vector describing the extrusion orientation
    used to find the local reference frame. As the sections are supposed to be
    placed in the XY plane, this vector will typically be (0,0,1).

    spine: 1D curve where the sections will be swept along. The spine-wise
    discretization is exactly the same as the node distribution of the spine.
    Hence, if you wish to control the discretization along the spine, change
    your spine preferrably with linelaw or splinelaw.

    sectionShapeLaw: It controls the law imposed to interpolate the airfoils
    between two sections. Can be ‘linear’, ‘nearest’, ‘zero’, ‘slinear’,
    ‘quadratic‘ or ‘cubic’.
    ______________
    Important note: We suppose that the input sections are positioned in the
    XY plane. The reference point from which the rotation
    angles are applied is the origin (0,0,0). If you wish to re-position
    this reference, you may apply a translation to the sections beforehand.
    '''

    # Checks the Number of sections and number of spine points
    Ns = len(sections)
    if Ns == 1:
        sections.append(sections[0])
        Ns += 1
    Ntot = spine[2]

    if not SpanPositions: SpanPositions = numpy.linspace(0,1,Ns)
    elif numpy.max(SpanPositions) > 1 or numpy.min(SpanPositions) < 0:
        print ('Aero sweepSections warning: SpanPositions badly imposed, \
            it shall be defined between [0,1]. Switching to uniform spacing.')
        SpanPositions = numpy.linspace(0,1,Ns)

    if not rotation:
        rotation = numpy.zeros(Ns)
    elif len(rotation) == 1:
        rotation = rotation[0]*numpy.ones(Ns)
    elif len(rotation) == 2:
        rotation = numpy.linspace(rotation[0],rotation[1],Ns)
    elif len(rotation) != Ns:
        raise ValueError("sweepSections: did not understand the imposed rotation: \
            number of sections: %d, number of rotation planes: %d"%(Ns,len(rotation)))


    # Invokes the Surface
    Surf = Converter.Array3D.convertArrays2Arrays3D(T.addkplane([sections[0]],Ntot-1))[0]

    # Makes the interpolation matrices based upon the provided sections
    W4dims = numpy.shape(Surf[1])
    InterpXmatrix = numpy.zeros((Ns,sections[0][2]),dtype=numpy.float64,order='F')
    InterpYmatrix = numpy.zeros((Ns,sections[0][2]),dtype=numpy.float64,order='F')
    for s in range(Ns):
        InterpXmatrix[s,:] = sections[s][1][0]
        InterpYmatrix[s,:] = sections[s][1][1]

    # Searches the elementary spanwise coordinate
    UnitCurvAbscissa = D.getDistribution(spine)[1][0]

    # Searches the tangent unitary vector of the spine
    Tang = D.getTangent(spine)

    # Makes the interpolation functions
    interpX = scipy.interpolate.interp1d( SpanPositions, InterpXmatrix, axis=0,
        kind=sectionShapeLaw, bounds_error=False, fill_value=InterpXmatrix[0,:])
    interpY = scipy.interpolate.interp1d( SpanPositions, InterpYmatrix, axis=0,
        kind=sectionShapeLaw, bounds_error=False, fill_value=InterpYmatrix[0,:])
    interpR = scipy.interpolate.interp1d( SpanPositions, rotation, axis=0,
        kind=rotationLaw, bounds_error=False, fill_value=rotation[0])

    # Constructs the surface plane by plane
    InterpolatedProfiles = [sections[0]]
    for k in range(Ntot):
        # X and Y coords are based upon the interpolating functions
        InterpolatedProfiles[k][1][0] = interpX(UnitCurvAbscissa[k])
        InterpolatedProfiles[k][1][1] = interpY(UnitCurvAbscissa[k])
        # Z coords results of a imposed twist + 3D rotation + translation
        # imposed twist
        InterpolatedProfiles[k] = T.rotate(InterpolatedProfiles[k],(0,0,0),
                                    (0.,0.,interpR(UnitCurvAbscissa[k])))
        # 3D rotation: We look for the local reference frame e1 e2 e3
        e3 = Tang[1][:,k].flatten()
        e2 = numpy.cross(e3,NormalDirection)
        e2 /= numpy.linalg.norm(e2) # normalization
        e1 = numpy.cross(e2,e3)
        InterpolatedProfiles[k] = T.rotate(InterpolatedProfiles[k],(0.,0.,0.),
                                                 ((1.,0.,0.),(0,1,0),(0,0,1)),(e1, e2, e3) )

        # translation
        InterpolatedProfiles[k] = T.translate(InterpolatedProfiles[k],
                                (spine[1][0][k],spine[1][1][k],spine[1][2][k]))
        # The interpolated profile at this point is done, we save it in our mesh
        # First we translate it into Arrays3D format
        InterpolatedProfiles[k] = Converter.Array3D.convertArrays2Arrays3D([InterpolatedProfiles[k]])[0]
        # and then we migrate the data
        Surf[1][0][:,0,k] = InterpolatedProfiles[k][1][0][:,:,0].flatten()
        Surf[1][1][:,0,k] = InterpolatedProfiles[k][1][1][:,:,0].flatten()
        Surf[1][2][:,0,k] = InterpolatedProfiles[k][1][2][:,:,0].flatten()
        # Invokes the next profile
        InterpolatedProfiles.append(sections[0])

    return Converter.Array3D.convertArrays3D2Arrays([Surf])[0]

