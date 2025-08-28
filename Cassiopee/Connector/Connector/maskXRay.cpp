/*    
    Copyright 2013 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

# include "connector.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Creation du masque X-Ray */
//=============================================================================
PyObject* K_CONNECTOR::maskXRay(PyObject* self, PyObject* args)
{
  PyObject* body;
  E_Float delta, tol;
  E_Int isNot, dim;
  if (!PYPARSETUPLE_(args, O_ R_ II_ R_,
                    &body, &delta, &dim, &isNot, &tol))
  {
      return NULL;
  }

  /* Extraction de la surface de masquage */
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;
  E_Int res = K_ARRAY::getFromArrays(
    body, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = unstrF.size();

  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "maskXRay: body arrays must be unstructured.");
    for (E_Int iu = 0; iu < nzones; iu++)
      RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
    return NULL;
  }

  if (nzones == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "maskXRay: at least one zone is required.");
    return NULL;
  }
  char* eltTypeb = NULL;

  // verification du type d'elements 
  for (E_Int i = 0; i < nzones; i++)
  {
    eltTypeb  = eltType[i];
    if (strcmp(eltTypeb, "TRI") != 0 && strcmp(eltTypeb, "BAR") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "maskXRay: body arrays must be all of TRI or BAR type.");
      for (E_Int iu = 0; iu < nzones; iu++)
        RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
      return NULL;
    }
  }
  
  // verification des coordonnees
  E_Int posxi, posyi, poszi;
  vector<E_Int> posxb; vector<E_Int> posyb; vector<E_Int> poszb;
  for (E_Int i = 0; i < nzones; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString[i]);
    if (posxi == -1 || posyi == -1 || poszi == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "maskXRay: body arrays must contain coordinates.");
      for (E_Int iu = 0; iu < nzones; iu++)
        RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
      return NULL;
    }
    posxi++; posyi++; poszi++;
    posxb.push_back(posxi); posyb.push_back(posyi); poszb.push_back(poszi);
  }

  E_Int dim1 = 1000; E_Int dim2 = 1000; 
  E_Int elevationDir = 3;
  if (strcmp(eltTypeb, "BAR") == 0 || dim == 2) elevationDir = 2;

  list<XRayPlane*> planes;
  E_Float xmin, ymin, zmin, xmax, ymax, zmax;
  E_Int ok = compCharacteristics(E_Int(isNot), elevationDir, dim1, dim2, 
                                 E_Float(tol), E_Float(delta),
                                 posxb, posyb, poszb, unstrF, cnt, planes,
                                 xmin, ymin, zmin, xmax, ymax, zmax);
  if (ok == -1) 
  {
    for (E_Int iu = 0; iu < nzones; iu++)
      RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
    PyErr_SetString(PyExc_TypeError,
                    "maskXRay: X-Ray mask has ambiguous faces.\n");
    return NULL;
  } 
  
  char varString[8]; strcpy(varString, "x,y,z");
  FldArrayI* connect = new FldArrayI();
  FldArrayF* coord = new FldArrayF();
  
  compXRayMaskInfo(elevationDir, planes, *coord);

  PyObject* tpl = K_ARRAY::buildArray(*coord, varString, *connect, 0);
  delete coord; delete connect;
  for (list<XRayPlane*>::iterator itr = planes.begin(); 
       itr != planes.end(); itr++)
  {
    delete [] (*itr)->tempZ;
    delete *itr;
  }
  planes.clear();
  for (E_Int iu = 0; iu < nzones; iu++)
    RELEASESHAREDU(obju[iu], unstrF[iu], cnt[iu]);
  return tpl;
}

//=========================================================================
/* Determination des points de percage du masque */
//=========================================================================
E_Int K_CONNECTOR::compCharacteristics(
  E_Int isNot, E_Int elevationDir, 
  E_Int dim1, E_Int dim2, 
  E_Float tol, E_Float delta,
  vector<E_Int>& posxt, 
  vector<E_Int>& posyt, 
  vector<E_Int>& poszt,
  vector<FldArrayF*>& fields, 
  vector<FldArrayI*>& cnt,
  list<XRayPlane*>& planes,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax)
{
  E_Int comp = 0;// nb de facettes ambigues
  // 1- calcul de la bb globale de fields
  K_COMPGEOM::globalBoundingBox(posxt, posyt, poszt, fields,
                                xmin, ymin, zmin, xmax, ymax, zmax);

  // 2- XRay plane structure
  XRayPlane* p = new struct XRayPlane;
  planes.push_back(p);
  
  list<XRayPlane*>::iterator itr = planes.begin();
  p = *itr;
  
  // Allocate each XRay planes
  p->xmin = xmin; p->ymin = ymin; p->zmin = zmin;
  p->xmax = xmax; p->ymax = ymax; p->zmax = zmax;
  if ( p->zmax < p->zmin+K_CONST::E_GEOM_CUTOFF ) p->zmax = p->zmin+K_CONST::E_GEOM_CUTOFF;
  p->ni = dim1; p->hi = (p->xmax - p->xmin)/(K_CONST::ONE*p->ni-1);
  if (elevationDir == 2) {p->nj = 1; p->hj = 1.;}
  else
  {p->nj = dim2;p->hj = (p->ymax - p->ymin)/(K_CONST::ONE*p->nj-1);}
  
  p->indir.malloc(p->ni * p->nj);
  p->tempZ = new vector<E_Float>[p->ni * p->nj];
  
  // Technique de perturbation des rayons pour resoudre les
  // facettes ambigues
  E_Int niter = 0;
  E_Int nitermax = 5;
  // Premier calcul
  FldArrayF epsilon(p->ni * p->nj, 2);
  epsilon.setAllValuesAtNull();
  
  // Compute Z intersections
  E_Int nzones = fields.size();
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    E_Int posxb = posxt[zone]; 
    E_Int posyb = posyt[zone]; 
    E_Int poszb = poszt[zone]; 
    FldArrayF& fieldb = *fields[zone];
    FldArrayI& cnb = *cnt[zone];
    comp = computeZ(elevationDir, xmin, ymin, xmax, ymax,
                    epsilon, fieldb.begin(posxb), fieldb.begin(posyb), 
                    fieldb.begin(poszb), cnb, p);
  }
  E_Int nambig = comp; comp = 0;
  while (nambig > 0 && niter <= nitermax)
  {
    // Compute Z intersections
    for (E_Int zone = 0; zone < nzones; zone++)
    {
      E_Int posxb = posxt[zone]; 
      E_Int posyb = posyt[zone]; 
      E_Int poszb = poszt[zone]; 
      FldArrayF& fieldb = *fields[zone];
      FldArrayI& cnb = *cnt[zone];
      comp = computeZ(elevationDir, xmin, ymin, xmax, ymax,
                      epsilon, fieldb.begin(posxb), fieldb.begin(posyb), 
                      fieldb.begin(poszb), cnb, p);
    }
    nambig = comp; comp = 0;
    niter++;
  }
  // Check if some faces are ambiguous
  comp = nambig;
  if (comp > 0) return -1;
  
  // Compact storage : reorder the Z array, enforce parity condition
  compactZ(isNot, elevationDir, delta, tol, 
           xmin, ymin, zmin, xmax, ymax, zmax, p);
  return 1;
}
//=============================================================================
/* Compute the Z pierce points for a body defined by fieldb */
//=============================================================================
E_Int K_CONNECTOR::computeZ(
  E_Int elevationDir,
  E_Float xmin, E_Float ymin, E_Float xmax, E_Float ymax, 
  FldArrayF& epsilon, 
  E_Float* xtb, E_Float* ytb, E_Float* ztb, 
  FldArrayI& cnb, struct XRayPlane* p)
{
  E_Int comp = 0;
  E_Int nelts = cnb.getSize();
  E_Int nvert = cnb.getNfld();
  E_Int* cn1 = cnb.begin(1);
  E_Int* cn2 = cnb.begin(2);  
  E_Int* cn3 = cnb.begin(1);

  if (nvert == 3) cn3 = cnb.begin(3);

//#pragma omp parallel for default(shared) reduction(+:comp)
  for (E_Int et = 0; et < nelts; et++)
  {
    E_Int ind1 = cn1[et]-1; E_Int ind2 = cn2[et]-1; E_Int ind3 = cn3[et]-1;
    comp = comp + triangleOnPlane(elevationDir, xmin, ymin, xmax, ymax, 
                                  epsilon, xtb, ytb, ztb, 
                                  ind1, ind2, ind3, p);
  }
  return comp;
}
//=============================================================================
E_Int K_CONNECTOR::triangleOnPlane(E_Int elevationDir, 
                                   E_Float xmina, E_Float ymina, 
                                   E_Float xmaxa, E_Float ymaxa,
                                   FldArrayF& epsilon, 
                                   E_Float* xt, E_Float* yt, E_Float* zt, 
                                   E_Int ind1, E_Int ind2, E_Int ind3,
                                   struct XRayPlane* p)
{
  E_Int comp = 0;
  E_Float eps = K_CONST::E_GEOM_CUTOFF;
  E_Float x, y, z;
  E_Float sgn1 = 0.;
  E_Float sgn2 = 0.;
  E_Int found;
  E_Float r1, r2, xp, yp;
  E_Int imin, imax, jmin, jmax;

  // Plane data
  E_Int ni = p->ni;
  E_Int nj = p->nj;
  E_Float xmin = p->xmin;
  E_Float ymin = p->ymin;
  E_Float hi = p->hi;
  E_Float hj = p->hj;
  E_LONG idum = -1;

  // Build triangle coordinates
  E_Float x0 = xt[ind1]; E_Float y0 = yt[ind1]; E_Float z0 = zt[ind1];
  E_Float x1 = xt[ind2]; E_Float y1 = yt[ind2]; E_Float z1 = zt[ind2];
  E_Float x2 = xt[ind3]; E_Float y2 = yt[ind3]; E_Float z2 = zt[ind3];
  // Find concerned rays
  imin = E_Int((x0-xmin-eps)/hi);
  imax = E_Int((x0-xmin+eps)/hi)+1;
  imin = K_FUNC::E_min( imin, E_Int((x1-xmin-eps)/hi) );
  imax = K_FUNC::E_max( imax, E_Int((x1-xmin+eps)/hi)+1 );
  imin = K_FUNC::E_min( imin, E_Int((x2-xmin-eps)/hi));
  imax = K_FUNC::E_max( imax, E_Int((x2-xmin+eps)/hi)+1);
  imax = K_FUNC::E_min( imax, ni-1 );
  imin = K_FUNC::E_max( imin, 0 );

  if (elevationDir != 2)
  {
    jmin = E_Int((y0-ymin-eps)/hj);
    jmax = E_Int((y0-ymin+eps)/hj)+1;
    jmin = K_FUNC::E_min( jmin, E_Int((y1-ymin-eps)/hj) );
    jmax = K_FUNC::E_max( jmax, E_Int((y1-ymin+eps)/hj)+1 );
    jmin = K_FUNC::E_min( jmin, E_Int((y2-ymin-eps)/hj) );
    jmax = K_FUNC::E_max( jmax, E_Int((y2-ymin+eps)/hj)+1 );
    jmax = K_FUNC::E_min( jmax, nj-1 );
    jmin = K_FUNC::E_max( jmin, 0 );
  }
  else
  {
    jmin = 0;
    jmax = 0;
  }

  for (E_Int j = jmin; j <= jmax; j++)
    for (E_Int i = imin; i <= imax; i++)
    {
      x = xmin + i * hi + epsilon(i+j*ni,1);
      y = ymin + j * hj + epsilon(i+j*ni,2);

      found = compIntersect(elevationDir,
                            x0, y0, z0, x1, y1, z1, x2, y2, z2,
                            x, y, p->zmin, p->zmax, z);

      switch (found)
      {
        case 1:
//#pragma omp critical
        {
          p->tempZ[i+j*ni].push_back(z);
        }
        break;
        case -1:
        case -2:
          sgn1 = K_NOISE::stdRand(&idum); sgn2 = K_NOISE::stdRand(&idum);
          if (sgn1 > 0.5) sgn1 = 0.1;
          else sgn1 = -0.1;
          if (sgn2 > 0.5) sgn2 = 0.1;
          else sgn2 = -0.1;
          r1 = K_NOISE::stdRand(&idum);
          r2 = K_NOISE::stdRand(&idum);
          xp = x + sgn1 * r1 * p->hi;
          if (xp < xmina) sgn1 = -sgn1;
          if (xp > xmaxa) sgn1 = -sgn1;
          yp = y + sgn2 * r2 * p->hj;
          if (yp < ymina) sgn2 = -sgn2;
          if (yp > ymaxa) sgn2 = -sgn2;
          // perturber les rayons d'un eps
          epsilon(i+j*ni,1) = sgn1 * r1 * p->hi;
          epsilon(i+j*ni,2) = sgn2 * r2 * p->hj;
          comp++;
          break;
        default:;
      }
      
    }
  return comp;
}

//-----------------------------------------------------------------------------
// Compute intersection between triangle defined by (x0, y0, z0), 
// (x1, y1, z1), (x2, y2, z2) and ray (x,y). Result in z if found.
// Return value is : 0 non-intersection
//                   1 intersection
//                  -1 ambiguous intersection, face is along z
//                  -2 intersection on a node of triangle
//-----------------------------------------------------------------------------
E_Int K_CONNECTOR::compIntersect(E_Int elevationDir,
                                 E_Float x0, E_Float y0, E_Float z0,
                                 E_Float x1, E_Float y1, E_Float z1,
                                 E_Float x2, E_Float y2, E_Float z2,
                                 E_Float x, E_Float y, 
                                 E_Float zmin, E_Float zmax,
                                 E_Float& z)
{
  E_Float delta, dx;
  E_Float r, t, l;
  E_Float alpha, beta;
  E_Float xmin, xmax;
  E_Float ymin, ymax;
  E_Float yy;//ordonnee du point d abscisse x sur la trace 
  E_Float delta1, delta2, delta3;
  E_Float eps = 1.e-10;
  E_Float eps2 = 1.e-7;//tolerance sur cos

  if (elevationDir == 3) // 3D case ---
  {
    E_Float x02 = x0-x2;
    E_Float y02 = y0-y2;
    E_Float z02 = z0-z2;
    E_Float x12 = x1-x2;
    E_Float y12 = y1-y2;
    E_Float z12 = z1-z2;
    E_Float dx2 = x-x2;
    E_Float dy2 = y-y2;
    E_Float dzmin2 = zmin-z2;
    E_Float dz = zmin-zmax;
    E_Float n012x = y02*z12 - z02*y12;
    E_Float n012y = z02*x12 - x02*z12;
    E_Float n012z = x02*y12 - y02*x12;
    
    delta = dz*n012z;
    delta1 = dz*(dx2*y12 - x12*dy2);
    delta2 = dz*(x02*dy2 - dx2*y02);
    delta3 = x02*(y12*dzmin2 - dy2*z12)
      + y02*(z12*dx2-dzmin2*x12)
      + z02*(x12*dy2-dx2*y12);

    E_Float n012  = 1./sqrt(n012x*n012x+n012y*n012y+n012z*n012z);    
    E_Float deltad = n012z * n012;
 
    if (K_FUNC::fEqualZero(deltad, eps2) != true)
    {
      E_Float deltai = 1./delta;
      r = delta1 * deltai;
      t = delta2 * deltai;
      l = delta3 * deltai;

      if (l <= 1.+eps && l >= -eps && r >= -eps && t >= -eps 
          && r+t <= 1.+2.*eps)
      {

        // ambiguite
        if (K_FUNC::fEqualZero(r, eps) == true || 
            K_FUNC::fEqualZero(t, eps) == true ||
            K_FUNC::fEqualZero(r+t-1., eps) == true)
        {
          return -1;
        }
        else
        {
          z = zmin - l * dz;
          return 1;
        }
      }
      else  return 0;
    }
    else // determinant = 0 
    {
      // det = 0-> trace de ABC tq les projections de AB et AC
      //sur le plan z=0  sont colineaires
      
      xmin = x0; ymin = y0;
      xmax = x0; ymax = y0;
      //test x1
      if (x1 < xmin)
      {
        xmin = x1; ymin = y1;
      }
      else 
      {
        xmax = x1; ymax = y1;
      }
      
      //test x2
      if (x2 < xmin)
      {
        xmin = x2; ymin = y2;
      }
      else 
      {
        xmax = x2; ymax = y2;
      }
      
      dx = xmax-xmin;
      // Les trois points  sont sur une droite (x,y)=cte
      if (K_FUNC::fEqualZero(dx, eps) == true)
      {
        if (K_FUNC::fEqualZero(x-xmin, eps) == true &&
            K_FUNC::fEqualZero(y-ymin, eps) == true)
          return -1;
        else return 0;
      }
      
      alpha = (ymax-ymin)/dx;
      beta  = ymin - alpha * xmin;

      // verifier que le rayon intersecte le segment trace de la facette
      if (x >= xmin-eps && x <= xmax+eps)
      {
        yy = alpha * x + beta;
        if (K_FUNC::fEqualZero(y-yy, eps) == true)
          return -1;
        else return 0;
      }
      return 0;
    }
  }
  else // 2D case, elevation in y ---
  {
    // in this case the Z elevation is in y
    // triangles Z coordinates are not meaningful
    delta = (x1 - x0);
    if (K_FUNC::fEqualZero(delta, eps) != true)
    {
      l = (x - x0)/delta;
    
      if (l <= 1.+eps && l >= -eps)
      {
        z = y0 + l*(y1 - y0);
        if (K_FUNC::fEqualZero(l, eps) == true ||
            K_FUNC::fEqualZero(l-1., eps) == true)
          return -2;
        else
          return 1;
      }
      else
      {
        return 0;
      }
    }
    else // determinant null
    {
      if (K_FUNC::fEqualZero(x - x0, eps) == true)
        return -1;
      else
        return 0;
    }
  } 
}

//=============================================================================
// Compact Z array
// Ensure that the number of pierce point on a ray is odd
// If not: add a pierce point next to the last one
// IN: tempZ
// OUT: indir + Z
//=============================================================================
void K_CONNECTOR::compactZ(
  E_Int isNot, E_Int elevationDir, E_Float delta, E_Float tol,
  E_Float xmin, E_Float ymin, E_Float zmin, 
  E_Float xmax, E_Float ymax, E_Float zmax,             
  XRayPlane* p)
{
  E_Int i, j, ind, in, start;
  E_Float z, temp;
  E_Int n = p->ni * p->nj;

  // Count the number of pierce points: size
  E_Int size = 0;
  for (i = 0; i < n; i++)
  {
    vector<E_Float>& v = p->tempZ[i];
    size = size+v.size();
  }

  E_Int s = 0;
  E_Int c = 0;
  p->Z.malloc(2*size); // dans le cas ou tous les pierce points sont impairs
  // Compacting
  c = 0;
  for (i = 0; i < n; i++)
  {
    // Compacting one ray
    start = c;
    p->indir[i] = start;
    vector<E_Float>& v = p->tempZ[i];
    s = v.size();
    for (j = 0; j < s; j++)
    {
      z = v[j];
      p->Z[c] = z;
      c++;
    }

    // Tri par ordre croissant
    for (j = start; j < start+s; j++)
    {
      // Find smallest element
      ind = j;
      for (in = j+1; in < start+s; in++)
      {
        if (p->Z[in] < p->Z[ind])
          ind = in;
      }
      if (ind != j)
      {
        temp = p->Z[j];
        p->Z[j] = p->Z[ind];
        p->Z[ind] = temp;
      }
    }

    // Elimination des doublons (points distant de moins de tol)
    c = start;
    for (j = start+1; j < start+s; j++)
    {
      if (K_FUNC::E_abs(p->Z[c] - p->Z[j]) > tol)
      {
        p->Z[c+1] = p->Z[j]; c++;
      }
    }
    if (s > 0) c++;
    // Ajout pour obtenir un nombre pair d'elements
    // Dans le cas classique, on ajoute un point au dessus du dernier point
    // d'intersection
    if ((c-start)%2 != 0)
    { 
      p->Z[c] = p->Z[c-1] + 1.e-6;
      c++;
    }
  }
  p->Z.reAlloc(c);

  // Add delta
  if (delta > 0.) addDeltaToZ(isNot, elevationDir, delta, 
                              xmin, ymin, zmin, xmax, ymax, zmax, p);

//   // Bilan
//   E_Int beg, end;
//   E_Int indp;
//   for (j = 0; j < p->nj; j++)
//     for (i = 0; i < p->ni; i++)
//     {
//       ind = i + p->ni * j;
//       indp = K_FUNC::E_min(ind+1, n-1);
//       beg = p->indir[ind];
//       end = p->indir[indp];
//       printf("Ray %f %d %d (%d) starts %d and ends %d.\n",
//       p->xmin+i*p->hi, i,j,ind, p->indir[ind], p->indir[indp]
//       for (E_Int l = beg; l < end; l++)
//         printf("%f ", p->Z[l]);
//       printf("\n");
//    }
}

//=============================================================================
// Add delta to intersection points ordinates
//=============================================================================
void K_CONNECTOR::addDeltaToZ(E_Int isNot, E_Int elevationDir, E_Float delta, 
                              E_Float xmin, E_Float ymin, E_Float zmin, 
                              E_Float xmax, E_Float ymax, E_Float zmax,       
                              XRayPlane* p)
{
  E_Int i, j, r;
  E_Int istart, iend;
  E_Int jstart;
  if (isNot == 1) delta = -delta; //masque inverse

  // Ajoute delta a la bounding box
  xmin = xmin - delta; ymin = ymin - delta;
  xmax = xmax + delta; ymax = ymax + delta;
  if (elevationDir != 2)
  {
    zmin = zmin - delta; zmax = zmax + delta;
  }
  E_Int n = p->ni * p->nj;

  istart = p->indir[0];
  jstart = p->indir[0];

  for (r = 0; r < n-1; r++)
  {
    iend = p->indir[r+1];

    j = jstart; i = istart;

    while (i < iend)
    {
      p->Z[j] = p->Z[i] - delta; i++; j++;
    
      while (i+1 < iend && p->Z[i] + delta >= p->Z[i+1] - delta)
        i = i+2;

      p->Z[j] = p->Z[i] + delta; j++; i++;
    }
    istart = iend;
    p->indir[r+1] = j;
    jstart = j;
  }

  // Last ray
  r = n-1;
  iend = p->Z.getSize();
  j = jstart; i = istart;

  while (i < iend)
  {
    p->Z[j] = p->Z[i] - delta; i++; j++;
    
    while (i+1 < iend && p->Z[i] + delta >= p->Z[i+1] - delta)
      i = i+2;

    p->Z[j] = p->Z[i] + delta; j++; i++;
  }
  p->Z.reAlloc(j);
}

//=============================================================================
/* mask info XRay */
//=============================================================================
void K_CONNECTOR::compXRayMaskInfo(E_Int elevationDir, 
                                   list<XRayPlane*>& planes, 
                                   FldArrayF& coord)
{
  list<XRayPlane*>::iterator itr;
  E_Int i, j, k;
  E_Int kstart, kend;
  E_Int size;
  E_Float x, y, z;
  E_Float hi, hj;

  E_Int nbTot = 0;
  for (itr = planes.begin(); itr != planes.end(); itr++)
  {
    XRayPlane& p = **itr;
    size = p.Z.getSize(); 
    nbTot = nbTot+size;
  }
  coord.malloc(nbTot,3);
  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);

  E_Int compt = 0;

  E_Int ni, nj;
    
  for (itr = planes.begin(); itr != planes.end(); itr++)
  {
    XRayPlane& p = **itr;
    hi = p.hi; hj = p.hj;
    ni = p.ni; nj = p.nj;
    x = p.xmin-hi; y = p.ymin-hj;

    if (elevationDir == 2)
    {
      for (i = 0 ; i < ni-1 ; i++)
      {
        x = x+hi;
        kstart = p.indir[i];
        kend = p.indir[i+1];
        for (k = kstart; k < kend; k++)
        {
          y = p.Z[k];
          z = 0.;
          xt[compt] = x;
          yt[compt] = y;
          zt[compt] = z;
          compt++;        
        } 
      }
      
      i = ni-1;
      x = x+hi;
      kstart = p.indir[i];
      kend = p.Z.getSize(); 
      for (k = kstart; k < kend; k++)
      {
        y = p.Z[k];
        z = 0.;
        xt[compt] = x;
        yt[compt] = y;
        zt[compt] = z;
        compt++;        
      } 
    }
    else // 3D
    {
      for (j = 0; j < nj; j++)
      {
        y = y+hj;
        x = p.xmin-hi;
        for (i = 0; i < ni; i++)
        {
          x = x+hi;
          kstart = p.indir[i+j*ni];
          
          if (i+j*ni == ni*nj-1)
            kend = p.Z.getSize();
          else 
            kend = p.indir[i+j*ni+1];

          for (k = kstart ; k < kend ; k++)
          {
            z = p.Z[k];
            xt[compt] = x;
            yt[compt] = y;
            zt[compt] = z;
            compt++;        
          }
        }
      }      
    }
  }
}
