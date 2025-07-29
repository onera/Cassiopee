import math, numpy

def GetReferencePointsLine(integrationDegree, quadratureType):

    # Gauss-Legendre Quadrature is a special case of Gauss-Jacobi for alpha=beta=0
    if quadratureType == "GaussLegendre":
        numPoints, points, weights = GaussJacobiGetPoints1D(integrationDegree, 0, 0);
    elif quadratureType == "GaussLobatto":
        numPoints,points,weights = GaussLobattoGetPoints1D(integrationDegree);
    # scale result from interval [-1,1] to [0,1]
    points = (points + 1.0) * 0.5;
    weights *= 0.5;
    return numPoints, points, weights;


def GaussJacobiGetPoints1D(integrationDegree, alpha, beta):

    # a rule with n points is exact for integrationDegree 2*n - 1
    #   we need to use the smallest integer n which yields at least the requested integrationDegree
    numPoints = (integrationDegree + 2) // 2;

    # compute points and weights in long double such that they are really accurate
    points = GaussJacobiComputeQuadraturePoints(numPoints, alpha, beta);
    weights = GaussJacobiComputeQuadratureWeights(numPoints, alpha, beta, points);

    return numPoints, points, weights;


def GaussJacobiComputeQuadraturePoints(numPoints, alpha, beta):
    points = numpy.zeros(numPoints);
    # select a tolerance, such that
    # a) it can be reached by Newton iteration
    # b) if possible, it is more accurate (two orders of magnitude) than FloatT relative precision
    tol =1e-16;
    maxIt = 100;

    # we determine the roots of the Jacobi polynomial of degree numPoints
    for i in range(numPoints):

        # initialize to roots of Chebyshev polynomial (alpha=-0.5, beta=-0.5), for which a closed expression exists
        xi = -math.cos((i + 0.5) * math.pi / numPoints);
        #// a better initial guess is the average of the Chebyshev point and the last point found (if it exists)
        if (i > 0):
            xi = 0.5 * (xi + points[i-1]);
        # Now apply Newton's method for root finding to the Polynomial P. In order to find a new root every time, we apply
        # polynomial deflation with the roots that have been found already, i.e. we search for the roots of
        # f = P / (\prod_j=0^{i-1}(x-x_j)), where f/f' = P/(P'-P*\sum_j=0^{i-1}(1/(x-x_j)))
        it = 0;
        for it in range(maxIt):
            sum_var = 0.0;
            for j in range(i):
                sum_var += 1.0 / (xi - points[j]);
            poly = GaussJacobiComputeJacobiPolynomial(numPoints, alpha, beta, xi);
            polyDeriv = GaussJacobiComputeJacobiPolynomialDerivative(numPoints, alpha, beta, xi);
            update = -poly / (polyDeriv - poly * sum_var);
            if (xi + update == xi):
                break;
            xi += update;
            if (abs(update) < tol and abs(poly) < tol):
                break;
        if (it == maxIt):
            raise ValueError("Failure to converge in determination of GaussJacobi quadrature points.");
        points[i] = xi;

    return points;

def GaussJacobiComputeQuadratureWeights(numPoints, alpha, beta, points):

    if (points.size != numPoints):
        raise ValueError("Number of provided points does not match requested number of weights");

    weights = numpy.zeros(numPoints);

    power = 2**(alpha + beta+1); # 2^(alpha+beta+1)
    factor = power * ComputeGammaQuotient(alpha + numPoints + 1, numPoints + 1) * ComputeGammaQuotient(beta + numPoints + 1, alpha + beta + numPoints + 1);

    for i in range(numPoints):

        x = points[i];
        deriv = GaussJacobiComputeJacobiPolynomialDerivative(numPoints, alpha, beta, x);
        weights[i] = factor / ((1.0 - x * x) * deriv * deriv);

    return weights;

def GaussJacobiComputeJacobiPolynomial(degree, alpha, beta, x):

    polyN = 1.0;

    if (degree == 0):
        return polyN;

    polyNPlusOne = 0.5 * ((alpha - beta) + (alpha + beta + 2) * x);

    # the definition is recursive, but we use a loop variant here
    for n in range(1,degree):
        polyNMinusOne = polyN;
        polyN = polyNPlusOne;

        tmp = 2 * n + alpha + beta;
        a1 = 2 * (n + 1) * (n + alpha + beta + 1) * tmp;
        a2 = (tmp + 1) * (alpha * alpha - beta * beta);
        a3 = tmp * (tmp + 1) * (tmp + 2);
        a4 = 2 * (n + alpha) * (n + beta) * (tmp + 2);

        polyNPlusOne = ((a2 + a3 * x) * polyN - a4 * polyNMinusOne) / a1;

    return polyNPlusOne;



def GaussJacobiComputeJacobiPolynomialDerivative(degree, alpha, beta, x):
    return 0.5 * (degree + alpha + beta + 1) * GaussJacobiComputeJacobiPolynomial(degree-1, alpha+1, beta+1, x);


def ComputeGammaQuotient(n, m):

    lower = max(1, min(n, m)-1);
    upper = max(n, m)-1;

    quotient = 1.0;
    for i in range(lower + 1, upper):
        quotient *= i;

    if (m > n):
        quotient = 1.0 / quotient;

    return quotient;





def GaussLobattoGetPoints1D(integrationDegree):

    # a rule with n points is exact for integration 2 * n - 3
    #   we need to use the smallest integer n which yields at least the requested integrationDegree
    numPoints = (integrationDegree + 4) // 2;
    points = numpy.zeros(numPoints);
    weights = numpy.zeros(numPoints);
    maxIt = 100;
    tol = 1e-16;

    points[0] = -1.0;
    points[numPoints - 1] = 1.0;
    weights[0] = 2.0 / ((numPoints - 1) * numPoints);
    weights[numPoints - 1] = weights[0];

    poly = numpy.zeros(numPoints);

    for i in range(1, numPoints//2):

        xi = -math.cos(i * math.pi / (numPoints - 1));
        it = 0;

        for it in range(maxIt):
            poly[0] = 1.0;
            poly[1] = xi;
            for idg in range(2,numPoints):
                poly[idg] = ((2 * idg - 1) * xi * poly[idg - 1] - (idg - 1) * poly[idg -2]) / idg;


            update = -(xi * poly[numPoints - 1] - poly[numPoints - 2]) / (numPoints * poly[numPoints - 1]);
            if (xi + update == xi):
                break;
            xi += update;
            if (abs(update) < tol * abs(xi)):
                break;

        if (it == maxIt):
            raise ValueError("Failure to converge in determination of GaussLobatto quadrature points.");

        points[i] = xi;
        weights[i] = (2.0 / ((numPoints - 1) * numPoints * poly[numPoints - 1] * poly[numPoints - 1]));
        points[numPoints - i - 1] = -points[i];
        weights[numPoints - i - 1] = weights[i];


    if ((numPoints + 1) % 2 == 0):
        xi = 0.0;
        poly[0] = 1.0;
        poly[1] = xi;
        for idg in range(2,numPoints):
            poly[idg] = ((2 * idg - 1) * xi * poly[idg - 1] - (idg - 1) * poly[idg -2]) / idg;

        points[(numPoints - 1) // 2] = 0.0;
        weights[(numPoints - 1) // 2] = (2.0 / ((numPoints - 1) * numPoints * poly[numPoints - 1] * poly[numPoints - 1]));
    return numPoints, points, weights;


def GaussLobattoGetPoints1D2(numPoints,rescale):
    integrationDegree = 2 * numPoints - 4;
    numPoints, points, weights = GaussLobattoGetPoints1D(integrationDegree);
    if rescale:
        points = (points + 1.0) * 0.5; # Rescale between 0 and 1
    return points;

def GetReferencePointsHexa(integrationDegree, quadratureType):

    # always use a general tensor product rule

    numPoints1D, points1D, weights1D = GetReferencePointsLine(integrationDegree, quadratureType);
    numPoints = numPoints1D * numPoints1D * numPoints1D;

    points = numpy.zeros((numPoints, 3));
    weights = numpy.zeros(numPoints);

    index = 0;
    for k in range(numPoints1D):
        for j in range(numPoints1D):
            for i in range(numPoints1D):
                points[index, 0] = points1D[i];
                points[index, 1] = points1D[j];
                points[index, 2] = points1D[k];
                weights[index] = weights1D[i] * weights1D[j] * weights1D[k];
                index += 1;

    assert(index == numPoints);
    return numPoints, points, weights;


def GetReferencePointsQuad(integrationDegree, quadratureType):

    # always use a general tensor product rule
    numPoints1D, points1D, weights1D = GetReferencePointsLine(integrationDegree, quadratureType);
    numPoints = numPoints1D * numPoints1D;

    points = numpy.zeros((numPoints, 2));
    weights = numpy.zeros(numPoints);

    index = 0;
    for j in range(numPoints1D):
        for i in range(numPoints1D):
            points[index, 0] = points1D[i];
            points[index, 1] = points1D[j];
            weights[index] = weights1D[i] * weights1D[j];
            index += 1;
    assert(index == numPoints);
    return numPoints, points, weights;

def GetReferencePointsData(integrationDegree, quadratureType, cellType):
    #points3 has 3 columns and variable number of rows (number of quadrature points per element)
    if cellType == 4: #(Quad)
        numPoints, points3, weigths = GetReferencePointsQuad(integrationDegree, quadratureType);
    elif cellType == 8: #Hexa
        numPoints, points3, weigths = GetReferencePointsHexa(integrationDegree, quadratureType);
    interpolationMatrix = GetInterpolationMatrix(cellType, points3);

    return weigths, interpolationMatrix

#@param cellData Array containing in the first three columns the physical coordinates of each element integration point.
#@param weights The weights of the reference points.
#@param nodalData The matrix consisting of the element's nodal data (first the coordinates, then velocities).
#@param[out] outElement Element metric data into which to put the integration points and weights.

def ComputeElementMetric(integrationDegree,quadratureType,cellType,nodalData):

    weights, interpolationMatrix = GetReferencePointsData(integrationDegree, quadratureType, cellType)
    numPoints = weights.size;
    cellData = (interpolationMatrix.dot(nodalData));
    return cellData


def GetInterpolationMatrix(cellType, referencePoints):
    if cellType == 4:
        numNodes = 4
    elif cellType == 8:
        numNodes = 8
    numPoints = referencePoints.shape[0];
    interpolationMatrix = numpy.zeros((numPoints, numNodes))

    for p in range(numPoints):
        coordRef = []
        coordRef.append(referencePoints[p, 0])
        coordRef.append(referencePoints[p, 1])
        if cellType==8:
            coordRef.append(referencePoints[p, 2])

        if cellType==4: shapeFunc = EvaluateShapeFunctionQuad(coordRef)
        elif cellType==8: shapeFunc = EvaluateShapeFunctionHexa(coordRef)
        for i in range(numNodes):
            interpolationMatrix[p, i] = shapeFunc[i]

    return interpolationMatrix

def EvaluateShapeFunctionHexa(coordRef):

    numNodes = 8
    permutation = [0, 1, 3, 2, 4, 5, 7, 6]

    order = 1 # order of the element;

    nNodes1D = order + 1;

    xiPoly   = numpy.zeros(nNodes1D);
    etaPoly  = numpy.zeros(nNodes1D);
    zetaPoly = numpy.zeros(nNodes1D);

    xiPoly   = ComputePolynomial1DComplete(coordRef[0], order);
    etaPoly  = ComputePolynomial1DComplete(coordRef[1], order);
    zetaPoly = ComputePolynomial1DComplete(coordRef[2], order);

    shapeFuncVal = numpy.zeros(numNodes);
    nodeindex = 0;

    for k in range(nNodes1D):
        for j in range(nNodes1D):
            for i in range(nNodes1D):
                shapeFuncVal[permutation[nodeindex]] = xiPoly[i] * etaPoly[j] * zetaPoly[k];
                nodeindex+=1;

    return shapeFuncVal

def EvaluateShapeFunctionQuad(coordRef):

    numNodes = 4 #nodes of the element, more if high order element
    permutation = [0, 1, 3, 2]

    order = 1 # order of the element;

    nNodes1D = order + 1;

    xiPoly   = numpy.zeros(nNodes1D);
    etaPoly  = numpy.zeros(nNodes1D);

    xiPoly   = ComputePolynomial1DComplete(coordRef[0], order);
    etaPoly  = ComputePolynomial1DComplete(coordRef[1], order);

    shapeFuncVal = numpy.zeros(numNodes);
    nodeindex = 0;

    for j in range(nNodes1D):
        for i in range(nNodes1D):
            shapeFuncVal[permutation[nodeindex]] = xiPoly[i] * etaPoly[j];
            nodeindex+=1;

    return shapeFuncVal


#/// Compute vector of the two-sided constructed one-dimensional Lagrangian interpolation polynomial contributions.
#/**
# * The contributions of the two-sided interpolation polynomials are used to compute the interpolation coefficients of
# * a complete one-dimensional line of support points with respect to a given one-dimensional interpolation point.
# * These contributions can be used to compute the interpolation coefficients of surface or volume cells by combining
# * the evaluations in each dimension of the element according to the element shape functions.
# * The two-sided constructed polynomials are used for directions of elements which are complete with respect to the
# * number of present support points, i.e. quads, hexahedra, the direction normal to the triangles of a prism.
# * The method may be used with an AD type to compute the corresponding derivatives.
# * @param x The one-dimensional node coordinate.
# * @param order Interpolation order, i.e. also order of the reference element.
# * @param[out] polynomials Pointer to polynomial contributions, must provide storage for at least order+1 T-values.
# */
def ComputePolynomial1DComplete(x, order):
    polynomials = numpy.zeros(order+1);
    if order==1:
        polynomials[0] = -x + 1.0;
        polynomials[1] =  x;
    elif order==2:
        xPow2 = x * x;
        polynomials[0] =  2.0 * xPow2 - 3.0 * x + 1.0;
        polynomials[1] = -4.0 * xPow2 + 4.0 * x;
        polynomials[2] =  2.0 * xPow2 - 1.0 * x;
    elif order==3:
        xPow2 = x * x;
        xPow3 = xPow2 * x;
        polynomials[0] =  -4.5 * xPow3 +  9.0 * xPow2 - 5.5 * x + 1.0;
        polynomials[1] =  13.5 * xPow3 - 22.5 * xPow2 + 9.0 * x;
        polynomials[2] = -13.5 * xPow3 + 18.0 * xPow2 - 4.5 * x;
        polynomials[3] =   4.5 * xPow3 -  4.5 * xPow2 + 1.0 * x;
    elif order==4:
        xPow2 = x * x;
        xPow3 = xPow2 * x;
        xPow4 = xPow3 * x;
        polynomials[0] =   (32.0/3.0) * xPow4 -  (80.0/3.0) * xPow3 +  (70.0/3.0) * xPow2 - (25.0/3.0) * x + 1.0;
        polynomials[1] = -(128.0/3.0) * xPow4 +       96.0  * xPow3 - (208.0/3.0) * xPow2 +      16.0  * x;
        polynomials[2] =        64.0  * xPow4 -      128.0  * xPow3 +       76.0  * xPow2 -      12.0  * x;
        polynomials[3] = -(128.0/3.0) * xPow4 + (224.0/3.0) * xPow3 - (112.0/3.0) * xPow2 + (16.0/3.0) * x;
        polynomials[4] =   (32.0/3.0) * xPow4 -       16.0  * xPow3 +  (22.0/3.0) * xPow2 -       1.0  * x;
    else:
        raise ValueError("Order not implemented!");
    return polynomials
