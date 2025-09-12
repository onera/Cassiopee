# test driver
import OCC.Driver as D

# Create parameter
hauteur = D.Scalar(1., name='hauteur')
hauteur.range = [0,1]

# Create points
P1 = D.Point( (0,0,0), name='P1' )

P2 = D.Point( (1,0,0), name='P2' )

P3 = D.Point( (1.5,1,0), name='P3' )
P3.y.range = [0,1]
D.Eq(P3.y.s, hauteur.s)

P4 = D.Point( (2.5,1,0), name='P4' )
P4.y.range = [0,1]
D.Eq(P4.y.s, hauteur.s)

P5 = D.Point( (3,0,0), name='P5' )

P6 = D.Point( (4,0,0), name='P6' )

P7 = D.Point( (4,2,0), name='P7' )

P8 = D.Point( (0,2,0), name='P8' )


# Create lines
spline1 = D.Spline1( [P1,P2,P3,P4,P5,P6], name='spline1' )

line1 = D.Line( P6, P7, name='line1' )
line2 = D.Line( P7, P8, name='line2' )
line3 = D.Line( P8, P1, name='line3' )

# Create sketch
sketch1 = D.Sketch([spline1, line1, line2, line3], name='sketch1')

# test
D.DRIVER.solve2()
D.DRIVER.instantiate([0.5])

sketch1.writeCAD('out.step')
