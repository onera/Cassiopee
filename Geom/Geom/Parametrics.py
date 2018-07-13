# A base of parametric curves and surfaces definition

base = {
    # curves (t parameter)
    'line': '{x} = {t} ; {y} = 0. ; {z} = 0.',
    'circle': '{x} = cos(2.*pi*{t}) ; {y} = sin(2.*pi*{t}) ;  {z} = 0.',
    # surfaces (t,u parameters )
    'plane': '{x} = {t} ; {y} = {u} ; {z} = 0.',
    'klein': '{x} = (3.*(1.+sin(2.*pi*{u})) + 2.*(1.-cos(2.*pi*{u})/2.)*cos(2.*pi*{t}))*cos(2.*pi*{u}) ; {y} = -2.*(1.-cos(2.*pi*{u})/2.)*sin(2.*pi*{t}) ; {z} = (4.+2.*(1.-cos(2.*pi*{u})/2.)*cos(2.*pi*{t}))*sin(2.*pi*{u})',
    'shell': '{x} = 1.2**(-pi/4+11*pi/4*{u})*(sin(pi*{t})**2 *sin(-pi/4+11*pi/4*{u})) ; {y} = 1.2**(-pi/4+11*pi/4*{u})*(sin(pi*{t})*cos(pi*{t})) ; {z} = 1.2**(-pi/4+11*pi/4*{u})*(sin(pi*{t})**2 *cos(-pi/4+11*pi/4*{u}))',
    'torus': '{x} = (1.+ 0.5*cos(2.*pi*{t}))*cos(2.*pi*{u}) ; {y} = 0.5*sin(2*pi*{t}) ; {z} = (1.+ 0.5*cos(2.*pi*{t}))*sin(2.*pi*{u})',
    'cosinus': '{x} = -1.+2.*{t} ; {y} = sin(pi*((-1.+2.*{t})**2+(-1+2.*{u})**2))/2. ; {z} = -1.+2.*{u}',
    'moebius': '{x} = cos(2.*pi*{u})+(-0.4+0.8*{t})*cos(2*pi*{u}/2)*cos(2*pi*{u}) ; {y} = (-0.4+0.8*{t})*sin(2*pi*{u}/2) ; {z} = sin(2*pi*{u})+(-0.4+0.8*{t})*cos(2*pi*{u}/2)*sin(2*pi*{u})',
    'riemann': '{x} = (-6.+12.*{t})*(-25.+50.*{u}) ; {y} = 30.*(-6.+12.*{t}) ; {z} = (-25.+50.*{u})**2 - (-6.+12.*{t})**2',
    'klein2': '{x} = (2. + cos(2*pi*{u}/2)* sin(2*pi*{t}) - sin(2*pi*{u}/2)* sin(2 *2*pi*{t}))* cos(2*pi*{u}) ; {y} = sin(2*pi*{u}/2)* sin(2*pi*{t}) + cos(2*pi*{u}/2) *sin(2* 2*pi*{t}) ; {z} = (2 + cos(2*pi*{u}/2)* sin(2*pi*{t}) - sin(2*pi*{u}/2)* sin(2 *2*pi*{t}))* sin(2*pi*{u})',
    'enneper': '{x} = (-2.+4.*{t}) -(-2.+4.*{t})**3/3.  + (-2+4*{t})*(-2+4*{u})**2; {y} = (-2+4*{t})**2 - (-2+4*{u})**2 ; {z} = (-2+4*{u}) -(-2+4*{u})**3/3  + (-2+4*{u})*(-2+4*{t})**2',
    'helix': '{x} = (1-0.1*cos(2*pi*{u}))*cos(4*pi*{t}) ; {y} = 0.1*(sin(2*pi*{u}) + 4*pi*{t}/1.7 -10) ; {z} = (1-0.1*cos(2*pi*{u}))*sin(4*pi*{t})',
    'hexahedron': '{x} = cos(2.*pi*{u})**3*cos(-1.3+2.6*{t})**3 ; {y} = sin(-1.3+2.6*{t})**3 ; {z} = sin(2.*pi*{u})**3*cos(-1.3+2.6*{t})**3',
    'sphere': '{x} = cos(-pi/2.+pi*{t})*cos(2.*pi*{u}) ; {y} = sin(-pi/2.+pi*{t}) ; {z} = cos(-pi/2.+pi*{t})*sin(2.*pi*{u})',
    'catalan': '{x} = (-pi+4.*pi*{t})-sin(-pi+4*pi*{t})*cosh(-2+4*{u}) ; {y} = 4*sin(1./2.*(-pi+4*pi*{t}))*sinh((-2+4*{u})/2) ; {z} = 1-cos(-pi+4*pi*{t})*cosh(-2+4*{u})',
    'toupie': '{x} = (abs(-1.+2.*{t})-1)**2 * cos(2*pi*{u}) ; {y} = -1.+2.*{t} ; {z} = (abs(-1.+2.*{t})-1)**2 * sin(2.*pi*{u})',
    'trumpet': '{x} = cos(2*pi*{t})*sin(0.03+1.97*{u}) ; {y} = (cos(0.03+1.97*{u})+log(tan(0.5*(0.03+1.97*{u})+0.00000001))) ; {z} = sin(2.*pi*{t})*sin(0.03+1.97*{u})',
    'bonbon': '{x} = 2*pi*{t} ; {y} = cos(2*pi*{t})*sin(2*pi*{u}) ; {z} = cos(2*pi*{t})*cos(2*pi*{u})',
    'helicoidal': '{x} = sinh(-pi+2.*pi*{u})*sin(-pi+2*pi*{t}) ; {y} = 3*(-pi+2*pi*{t}) ; {z} = -sinh(-pi+2*pi*{u})*cos(-pi+2*pi*{t})',
    'horn': '{x} = (2 + {t}*cos(2.*pi*{u}))*sin(2*pi*{t}) ; {y} = {t} *sin(2*pi*{u}) ; {z} = (2 + {t}*cos(2*pi*{u}))*cos(2*pi*{t}) + 2*{t}',
    'steiner': '{x} = sin(2.*(-pi/2.+pi*{t})) * cos(-pi/2.+pi*{u}) * cos(-pi/2.+pi*{u}) ; {y} = cos(-pi/2.+pi*{t}) * sin(2*(-pi/2.+pi*{u})) ; {z} = sin(-pi/2.+pi*{t}) * sin(2.*(-pi/2.+pi*{u}))',
    'crosscap': '{x} = sin(-pi/2.+pi*{t}) * sin(-pi+2*pi*{u}) / 2 ; {y} = cos(-pi+2*pi*{t}) * cos(-pi/2+pi*{u}) * cos(-pi/2.+pi*{u}) ; {z} = sin(-pi+2.*pi*{t}) * cos(-pi/2.+pi*{u}) * cos(-pi/2.+pi*{u})',
    'cliffordtorus': '{x} = cos(pi*{t}+2.*pi*{u})/(sqrt(2.)+cos(2.*pi*{u}-pi*{t})) ; {y} = sin(2.*pi*{u}-pi*{t})/(sqrt(2.)+cos(2*pi*{u}-pi*{t})) ; {z} = sin(pi*{t}+2*pi*{u})/(sqrt(2.)+cos(2.*pi*{u}-pi*{t}))',
    'enneper2': '{x} = 1.2*{t}*cos(-pi+2.*pi*{u})-1.2*{t}**3/3.*cos(3*(-pi+2.*pi*{u})) ; {y} = (1.2*{t})**2*cos(2.*(-pi+2*pi*{u})) ; {z} = -1.2*{t}*sin(-pi+2*pi*{u})-(1.2*{t})**(3)/3.*sin(3*(-pi+2.*pi*{u}))',
    'dini': '{x} = cos(12.4*{t})*sin(0.1+1.9*{u}) ; {y} = (cos(0.1+1.9*{u})+log(tan((0.1+1.9*{u})/2))) + 0.2*12.4*{t} ; {z} = sin(12.4*{t})*sin(0.1+1.9*{u})',
    'drop': '{x} = 2.*{t}*cos(2.*pi*{u}) ; {y} = exp(-4.*{t}*{t})*(sin(4.*pi*{t}) - 2.*{t}*cos(6.*pi*{u})) ; {z} = 2.*{t}*sin(2.*pi*{u})'
    }

