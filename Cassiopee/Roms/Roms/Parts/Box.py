# Box part
import Roms.Driver as D

def createPart(name):

    T1 = D.Part(name)

    # Create parameters
    height = T1.Scalar('height')
    height.range = [0., 2., 0.1]

    width = T1.Scalar('width')
    width.range = [0., 2., 0.1]

    depth = T1.Scalar('depth')
    depth.range = [0., 2., 0.1]

    # Create points
    P1 = T1.Point('P1', (0,0,0))
    T1.Eq(P1.x, 0.)
    T1.Eq(P1.y, 0.)
    T1.Eq(P1.z, 0.)

    P2 = T1.Point('P2', (1,0,0))
    T1.Eq(P2.x, P1.x + width)

    P3 = T1.Point('P3', (1,1,0))
    T1.Eq(P3.x, P1.x + width)
    T1.Eq(P3.y, P1.y + height)
    
    P4 = T1.Point('P4', (0,1,0))
    T1.Eq(P4.y, P1.y + height)

    P5 = T1.Point('P5', (0,0,1))
    T1.Eq(P5.z, P1.z + depth)

    P6 = T1.Point('P6', (1,0,1))
    T1.Eq(P6.x, P1.x + width)
    T1.Eq(P6.z, P1.z + depth)

    P7 = T1.Point('P7', (1,1,1))
    T1.Eq(P7.x, P1.x + width)
    T1.Eq(P7.y, P1.y + height)
    T1.Eq(P7.z, P1.z + depth)
    
    P8 = T1.Point('P8', (0,1,1))
    T1.Eq(P8.y, P1.y + height)
    T1.Eq(P8.z, P1.z + depth)
    
    # Create lines
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    line5 = T1.Line('line5', P5, P6)
    line6 = T1.Line('line6', P6, P7)
    line7 = T1.Line('line7', P7, P8)
    line8 = T1.Line('line8', P8, P5)
    line9 = T1.Line('line9', P1, P5)
    line10 = T1.Line('line10', P2, P6)
    line11 = T1.Line('line11', P3, P7)
    line12 = T1.Line('line12', P4, P8)

    # Create sketch
    sketch1 = T1.Sketch('sketch1',
                        [line1, line2, line3, line4],
                        h=[0.01,0.01,0.01])
    face1 = T1.Fill('face1', sketch1, reverse=1)

    sketch2 = T1.Sketch('sketch2',
                        [line5, line6, line7, line8],
                        h=[0.01,0.01,0.01])
    face2 = T1.Fill('face2', sketch2)

    sketch3 = T1.Sketch('sketch3',
                        [line1, line10, line5, line9],
                        h=[0.01,0.01,0.01])
    face3 = T1.Fill('face3', sketch3)

    sketch4 = T1.Sketch('sketch4',
                        [line2, line11, line6, line10],
                        h=[0.01,0.01,0.01])
    face4 = T1.Fill('face4', sketch4)

    sketch5 = T1.Sketch('sketch5',
                        [line3, line11, line7, line12],
                        h=[0.01,0.01,0.01])
    face5 = T1.Fill('face5', sketch5)

    sketch6 = T1.Sketch('sketch6',
                        [line4, line12, line8, line9],
                        h=[0.01,0.01,0.01])
    face6 = T1.Fill('face6', sketch6)

    surface = T1.Merge('surface', [face1,face2,face3,face4,face5,face6])
    
    # solve result
    T1.solve()

    # example of instantiation
    #T1.instantiate({'width':1., 'height':1., 'depth':1.})

    return T1