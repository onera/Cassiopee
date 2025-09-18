/*    
    Copyright 2013-2025 Onera.

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

# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

#define NORM2(vx,vy,vz) vx*vx+vy*vy+vz*vz
#define SCAL(v1x,v1y,v1z,v2x,v2y,v2z) v1x*v2x+v1y*v2y+v1z*v2z

// ============================================================================
/* Smooth a mesh */
// ============================================================================
PyObject* K_TRANSFORM::smooth(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float eps;
  E_Int niter, type;
  PyObject* fixedConstraint; PyObject* projConstraint;
  E_Float delta; E_Float xR, yR, zR; E_Float radius;

  if (!PYPARSETUPLE_(args, O_ R_ II_ OO_ R_ TRRR_ R_,
                    &array, &eps, &niter, &type,
                    &fixedConstraint, &projConstraint, &delta,
                    &xR, &yR, &zR, &radius))
  {
    return NULL;
  }
  // Check array
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 =  K_ARRAY::getFromArray3(array, varString1, 
                                       f1, im1, jm1, km1, cn1, eltType1);

  if (res1 != 2)
  {
    if (res1 == 1) RELEASESHAREDS(array, f1);
    PyErr_SetString(PyExc_TypeError,
                    "smooth: array must be unstructured.");
    return NULL;
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
   
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDU(array, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "smooth: can't find coordinates in array.");
    return NULL;
  }
  posx1++; posy1++; posz1++;

  // fixed constraint
  E_Int fixedConstraintOn = 1;
  if (PyList_Check(fixedConstraint) == true && PyList_Size(fixedConstraint) == 0)
    fixedConstraintOn = 0;

  E_Int im2, jm2, km2;
  FldArrayF* f2=NULL; FldArrayI* cn2=NULL;
  char* varString2; char* eltType2; 
  E_Int res2 = 0;
  if (fixedConstraintOn == 1)
  {
    res2 = K_ARRAY::getFromArray3(fixedConstraint, varString2, 
                                  f2, im2, jm2, km2, cn2, eltType2);
    
  }
  if (fixedConstraintOn == 1 && res2 != 2)
  {
    RELEASESHAREDU(array, f1, cn1);
    if (res2 == 1) RELEASESHAREDS(fixedConstraint, f2);
    PyErr_SetString(PyExc_TypeError,
                    "smooth: constraints must be unstructured.");
    return NULL;
  }
  
  E_Int posx2=0, posy2=0, posz2=0;
  if (fixedConstraintOn == 1)
  {
    posx2 = K_ARRAY::isCoordinateXPresent(varString2);
    posy2 = K_ARRAY::isCoordinateYPresent(varString2);
    posz2 = K_ARRAY::isCoordinateZPresent(varString2);
   
    if (posx2 == -1 || posy2 == -1 || posz2 == -1)
    {
      RELEASESHAREDU(array, f1, cn1); RELEASESHAREDU(fixedConstraint, f2, cn2);
      PyErr_SetString(PyExc_TypeError,
                      "smooth: can't find coordinates in constraint array.");
      return NULL;
    }
    posx2++; posy2++; posz2++;
  }

  // proj constraint
  E_Int projConstraintOn = 1;
  if (PyList_Check(projConstraint) == true && PyList_Size(projConstraint) == 0)
    projConstraintOn = 0;

  E_Int im3, jm3, km3;
  FldArrayF* f3=NULL; FldArrayI* cn3=NULL;
  char* varString3; char* eltType3;
  E_Int res3 = 0;
  if (projConstraintOn == 1)
  {
    res3 = K_ARRAY::getFromArray3(projConstraint, varString3, 
                                  f3, im3, jm3, km3, cn3, eltType3);
  }
  if (projConstraintOn == 1 && res3 != 2)
  {
    RELEASESHAREDU(array, f1, cn1); 
    if (fixedConstraintOn == true) RELEASESHAREDU(fixedConstraint, f2, cn2);
    if (res3 == 1) RELEASESHAREDS(projConstraint, f3);
    PyErr_SetString(PyExc_TypeError,
                    "smooth: constraints must be unstructured.");
    return NULL;
  }

  E_Int posx3=0, posy3=0, posz3=0;
  if (projConstraintOn == 1)
  {
    posx3 = K_ARRAY::isCoordinateXPresent(varString3);
    posy3 = K_ARRAY::isCoordinateYPresent(varString3);
    posz3 = K_ARRAY::isCoordinateZPresent(varString3);
   
    if (posx3 == -1 || posy3 == -1 || posz3 == -1)
    {
      RELEASESHAREDU(array, f1, cn1); 
      if (fixedConstraintOn == true) RELEASESHAREDU(fixedConstraint, f2, cn2);
      RELEASESHAREDU(projConstraint, f3, cn3);
      PyErr_SetString(PyExc_TypeError,
                      "smooth: can't find coordinates in constraint array.");
      return NULL;
    }
    posx3++; posy3++; posz3++;
  }

  // array
  E_Int csize;
  if (strcmp(eltType1, "NGON") == 0) csize = cn1->getSize()*cn1->getNfld();
  else csize = 1;
  PyObject*  tpl = K_ARRAY::buildArray(
    f1->getNfld(), varString1, 
    f1->getSize(), cn1->getSize(), -1, eltType1, false, csize);

  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cn1->begin(), cn1->getSize()*cn1->getNfld());
  FldArrayI cnn(cn1->getSize(), cn1->getNfld(), cnnp, true);
  E_Float* coordop = K_ARRAY::getFieldPtr(tpl);
  memcpy(coordop, f1->begin(), f1->getSize()*f1->getNfld()*sizeof(E_Float));
  FldArrayF coordo(f1->getSize(), f1->getNfld(), coordop, true);

  E_Int npts = f1->getSize();
  // Pour le stockage iteratif
  FldArrayF* coord = new FldArrayF(npts, 3);
  coord->setOneField(*f1, posx1, 1);
  coord->setOneField(*f1, posy1, 2);
  coord->setOneField(*f1, posz1, 3);
  // Pour stocker la projection sur la contrainte
  FldArrayF* proj = new FldArrayF(npts, 3);
  E_Float* projx=NULL, *projy=NULL, *projz=NULL;
  // Tableau indiquant les pts que l'on peut bouger
  FldArrayF* move = new FldArrayF(npts);
  E_Float* movep = move->begin();
  move->setAllValuesAt(1.);

  // Lissage isotrope / umbrella avec les vertex
  //E_Int nit = 0;
  vector< vector<E_Int> > cVN(npts);// vertex/elt
  if (strcmp(eltType1, "NGON") == 0) K_CONNECT::connectNG2VNbrs(*cn1, cVN);
  else K_CONNECT::connectEV2VNbrs(*cn1, cVN);    

  E_Float* cx = coord->begin(1);
  E_Float* cy = coord->begin(2);
  E_Float* cz = coord->begin(3);
  E_Float* f1x = coordo.begin(posx1);
  E_Float* f1y = coordo.begin(posy1);
  E_Float* f1z = coordo.begin(posz1);
   
  // Les points proches de la contrainte ne peuvent bouger
  if (fixedConstraintOn == 1)
  {
    proj->setOneField(coordo, posx1, 1);
    proj->setOneField(coordo, posy1, 2);
    proj->setOneField(coordo, posz1, 3);
    projx = proj->begin(1);
    projy = proj->begin(2);
    projz = proj->begin(3);
      
    // Projection du maillage sur les contraintes fixes (seulement si
    // la surface est TRI ou BAR)
    if (cn2->getNfld() == 2 || cn2->getNfld() == 3)
    {
      K_COMPGEOM::projectOrthoWithPrecond(npts, *cn2,
                                          posx2, posy2, posz2,
                                          *f2, projx, projy, projz);
    }
    else
    {
      // Si on ne peut pas projeter, on cherche par le point
      // le plus proche (par KdTree)
      ArrayAccessor<FldArrayF> coordAcc(*f2, posx2, posy2, posz2);
      KdTree<FldArrayF> kdt(coordAcc, E_EPSILON);
      E_Float* f2x = f2->begin(posx2);
      E_Float* f2y = f2->begin(posy2);
      E_Float* f2z = f2->begin(posz2);
      E_Float pt[3];
      for (E_Int ind = 0; ind < npts; ind++)
      {
        E_Int indw2;
        pt[0] = f1x[ind]; pt[1] = f1y[ind]; pt[2] = f1z[ind];
        indw2 = kdt.getClosest(pt);
        projx[ind] = f2x[indw2];
        projy[ind] = f2y[indw2];
        projz[ind] = f2z[indw2];
      }
    }

#pragma omp parallel default(shared)
    {
      E_Float dx, dy, dz, distv;
      E_Int nbV; 
      E_Float deltar; E_Int nov; E_Float rayon;
      
#pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        vector<E_Int>& v = cVN[ind]; // vertex voisins
        nbV = v.size();
        rayon = 0.;
        for (E_Int vi = 0; vi < nbV; vi++)
        {
          nov = v[vi]-1;
          dx = cx[ind]-cx[nov];
          dy = cy[ind]-cy[nov];
          dz = cz[ind]-cz[nov];
          //rayon = K_FUNC::E_max(dx*dx+dy*dy+dz*dz, rayon);
          rayon += dx*dx+dy*dy+dz*dz;
        }
        rayon = rayon / nbV;
     
        dx = projx[ind]-cx[ind];
        dy = projy[ind]-cy[ind];
        dz = projz[ind]-cz[ind];
        distv = dx*dx + dy*dy + dz*dz;
        deltar = delta*rayon;
        if (distv < deltar) {
          // contrainte lineaire
          movep[ind] = distv/deltar;
          // contrainte quadratique
          //movep[ind] = -(distv*distv)/(deltar*deltar)+2*distv/deltar;
        }
      }
    }
  }
  
  if (type == 0) // isotrope umbrella -> force la regularite des volumes
  {
    umbrella(*coord, coordo, *move,
             *proj, f3, cn3,
             posx1, posy1, posz1,
             posx3, posy3, posz3,
             projConstraintOn, delta, type, eps, niter, cVN,
             xR, yR, zR, radius);
  }
  else if (type == 1) // scaled umbrella
  {
    umbrella(*coord, coordo, *move,
             *proj, f3, cn3,
             posx1, posy1, posz1,
             posx3, posy3, posz3,
             projConstraintOn, delta, type, eps/5., niter, cVN,
             xR, yR, zR, radius);
  }
  else if (type == 2) // taubin lambda/mu
  {
    for (E_Int nit = 0; nit < niter; nit++)
    {
      umbrella(*coord, coordo, *move,
               *proj, f3, cn3,
               posx1, posy1, posz1,
               posx3, posy3, posz3,
               projConstraintOn, delta, 0, -0.333, 1, cVN,
               xR, yR, zR, radius);
      umbrella(*coord, coordo, *move,
               *proj, f3, cn3,
               posx1, posy1, posz1,
               posx3, posy3, posz3,
               projConstraintOn, delta, 0, +0.333, 1, cVN,
               xR, yR, zR, radius);
    }
  }

  // Build array
  delete coord; delete proj; delete move;
  RELEASESHAREDU(array, f1, cn1); 
  if (fixedConstraintOn == true) RELEASESHAREDU(fixedConstraint, f2, cn2);
  if (projConstraintOn == true) RELEASESHAREDU(projConstraint, f3, cn3);
  return tpl;
}

//==============================================================================
// IN/OUT: coord
// IN: coordo: stockage temporaire
// IN: move: si move[ind]=0, le pt ne bouge pas
// IN: proj: stockage temporaire du projete
// IN: f3: contrainte de glissement (projection)
// IN: cn3: connectivite de la contrainte de glissement
// IN: posx1, posy1, posz1: positions des coordonnees dans coord
// IN: posx3, posy3, posz3: positions des coordonnees dans f3
// IN: projConstraintOn: =1 (on projete sur la contrainte)
// IN: delta: force de la contrainte
// IN: type: =0 (isotrope), =1 (scale)
// IN: eps: coefficient d'integration
// IN: niter: nbre d'iteration de lissage
// IN: cVN: connectivite vertex-vertex voisins
// IN: xR,yR,zR: si local smoothing = center
// IN: radius: si local smoothing = radius
//==============================================================================
void  K_TRANSFORM::umbrella(FldArrayF& coord, FldArrayF& coordo,
                            FldArrayF& move, FldArrayF& proj,
                            FldArrayF* f3, FldArrayI* cn3,
                            E_Int posx1, E_Int posy1, E_Int posz1,
                            E_Int posx3, E_Int posy3, E_Int posz3,
                            E_Int projConstraintOn, E_Float delta, E_Int type,
                            E_Float eps, E_Int niter,
                            vector< vector<E_Int> > &cVN,
                            E_Float xR, E_Float yR, E_Float zR, E_Float radius)
{
  E_Int nit = 0;
  E_Float* cx = coord.begin(1);
  E_Float* cy = coord.begin(2);
  E_Float* cz = coord.begin(3);
  E_Float* f1x = coordo.begin(posx1);
  E_Float* f1y = coordo.begin(posy1);
  E_Float* f1z = coordo.begin(posz1);
  E_Float* movep = move.begin();
  E_Int npts = coord.getSize();

  while (nit < niter)
  {
#pragma omp parallel default(shared)
    {
      E_Float* projx=NULL, *projy=NULL, *projz=NULL;
      E_Float /*alpha,*/ eps2;
      E_Float deltax, deltay, deltaz, dx, dy, dz, dist, /*distv,*/ w, sum;
      E_Int nbV, nov;
      E_Float loc, r; 

      if (type == 0 && radius < 0) // isotrope
      {
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          vector<E_Int>& v = cVN[ind]; // vertex voisins
          nbV = v.size();
          dx = 0.; dy = 0.; dz = 0.;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            nov = v[vi]-1;
            dx += cx[nov]-cx[ind];
            dy += cy[nov]-cy[ind];
            dz += cz[nov]-cz[ind];
          }
          sum = nbV;
          eps2 = movep[ind]*eps/sum;
          f1x[ind] = cx[ind] + eps2 * dx;
          f1y[ind] = cy[ind] + eps2 * dy;
          f1z[ind] = cz[ind] + eps2 * dz;
        }
      }
      else if (type == 1 && radius < 0) // scale
      {
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          vector<E_Int>& v = cVN[ind]; // vertex voisins
          nbV = v.size();
          deltax = 0.; deltay = 0.; deltaz = 0.; sum = 0.;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            nov = v[vi]-1;
            dx = cx[nov]-cx[ind];
            dy = cy[nov]-cy[ind];
            dz = cz[nov]-cz[ind];
            w = sqrt(dx*dx+dy*dy+dz*dz);
            if (w > 1.e-10) w = 1./w;
            else w = 1.e10;
            deltax += w*dx;
            deltay += w*dy;
            deltaz += w*dz;
            sum += w;
          }
          eps2 = movep[ind]*eps/sum;
          f1x[ind] = cx[ind] + eps2 * deltax;
          f1y[ind] = cy[ind] + eps2 * deltay;
          f1z[ind] = cz[ind] + eps2 * deltaz;
        }
      }
      else if (type == 0 && radius >= 0) // isotrope
      {
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          dx = cx[ind]-xR;
          dy = cy[ind]-yR;
          dz = cz[ind]-zR;
          r = (dx*dx+dy*dy+dz*dz) / (radius*radius); 
          loc = exp(-r);

          vector<E_Int>& v = cVN[ind]; // vertex voisins
          nbV = v.size();
          dx = 0.; dy = 0.; dz = 0.;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            nov = v[vi]-1;
            dx += cx[nov]-cx[ind];
            dy += cy[nov]-cy[ind];
            dz += cz[nov]-cz[ind];
          }
          sum = nbV;
          eps2 = loc*movep[ind]*eps/sum;
          f1x[ind] = cx[ind] + eps2 * dx;
          f1y[ind] = cy[ind] + eps2 * dy;
          f1z[ind] = cz[ind] + eps2 * dz;
        }
      }
      else if (type == 1 && radius >= 0) // scale
      {
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          dx = cx[ind]-xR;
          dy = cy[ind]-yR;
          dz = cz[ind]-zR;
          r = (dx*dx+dy*dy+dz*dz) / (radius*radius); 
          loc = exp(-r);

          vector<E_Int>& v = cVN[ind]; // vertex voisins
          nbV = v.size();
          deltax = 0.; deltay = 0.; deltaz = 0.; sum = 0.;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            nov = v[vi]-1;
            dx = cx[nov]-cx[ind];
            dy = cy[nov]-cy[ind];
            dz = cz[nov]-cz[ind];
            w = sqrt(dx*dx+dy*dy+dz*dz);
            if (w > 1.e-10) w = 1./w;
            else w = 1.e10;
            deltax += w*dx;
            deltay += w*dy;
            deltaz += w*dz;
            sum += w;
          }
          eps2 = loc*movep[ind]*eps/sum;
          f1x[ind] = cx[ind] + eps2 * deltax;
          f1y[ind] = cy[ind] + eps2 * deltay;
          f1z[ind] = cz[ind] + eps2 * deltaz;
        }
      }

      if (projConstraintOn == 1)
      {
        proj.setOneField(coordo, posx1, 1);
        proj.setOneField(coordo, posy1, 2);
        proj.setOneField(coordo, posz1, 3);
        projx = proj.begin(1);
        projy = proj.begin(2);
        projz = proj.begin(3);
      
        // Projection du maillage sur les contraintes
        K_COMPGEOM::projectOrthoWithPrecond(npts, *cn3,
                                            posx3, posy3, posz3,
                                            *f3, projx, projy, projz);

#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          vector<E_Int>& v = cVN[ind]; // vertex voisins
          E_Float nbV = v.size();
          dx = projx[ind]-cx[ind];
          dy = projy[ind]-cy[ind];
          dz = projz[ind]-cz[ind];
          E_Float distP = sqrt(dx*dx + dy*dy + dz*dz);
          distP = K_FUNC::E_max(distP, 1.e-10);
          //E_Bool edge = false;

          /* Cherche le rayon le plus proche de la direction de proj */
          E_Float v1x = dx/distP; E_Float v1y = dy/distP; E_Float v1z = dz/distP;
          E_Float maxScal = -1.e6; E_Float maxL = -1.;
          E_Float v2x, v2y, v2z, s;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            E_Int nov = v[vi]-1;
            dx = cx[nov]-cx[ind];
            dy = cy[nov]-cy[ind];
            dz = cz[nov]-cz[ind];
            dist = sqrt(dx*dx + dy*dy + dz*dz);
            dist = K_FUNC::E_max(dist, 1.e-10);
            v2x = dx/dist; v2y = dy/dist; v2z = dz/dist;
            s = SCAL(v1x,v1y,v1z,v2x,v2y,v2z);
            if (s > maxScal) { maxScal = s; maxL = dist; }
          }
          //printf("%d: %f %f\n",ind,distP,maxL);
          //printf("delta %f\n", delta);

          // verifie si les vertex voisins ont deja ete projetes
          E_Int dejaProj = 0;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            E_Int nov = v[vi]-1;
            if ((f1x[nov]-projx[nov])*(f1x[nov]-projx[nov])+
                (f1y[nov]-projy[nov])*(f1y[nov]-projy[nov])+
                (f1z[nov]-projz[nov])*(f1z[nov]-projz[nov]) < 0.1*maxL*maxL) dejaProj++;
          }

          if (distP < 0.25*maxL) // on projete exactement
          {
            //printf("projection %f %f\n", distP, maxL);
            f1x[ind] = projx[ind]; f1y[ind] = projy[ind]; f1z[ind] = projz[ind];
          }
          else if (distP < 0.8*maxL) // on limite la projection
          {
            //printf("limited projection %f %f\n", distP, maxL);
            // Essai pour projeter plus les pts proches
            //f1x[ind] = f1x[ind]+eps*(0.8*maxL-distP)*v1x;
            //f1y[ind] = f1y[ind]+eps*(0.8*maxL-distP)*v1y;
            //f1z[ind] = f1z[ind]+eps*(0.8*maxL-distP)*v1z;
            // Reglage pas mal pour eviter trop de proj
            f1x[ind] = f1x[ind]+0.4*eps*maxL*v1x;
            f1y[ind] = f1y[ind]+0.4*eps*maxL*v1y;
            f1z[ind] = f1z[ind]+0.4*eps*maxL*v1z;
          }

          /*
          // Calcul le rayon de la cellule - moyen ou min?
          //E_Float rayon = 0.;
          E_Float rayon = 1.e6;
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            E_Int nov = v[vi]-1;
            dx = cx[ind]-cx[nov];
            dy = cy[ind]-cy[nov];
            dz = cz[ind]-cz[nov];
            rayon = K_FUNC::E_max(dx*dx+dy*dy+dz*dz, rayon);
            rayon += dx*dx+dy*dy+dz*dz;
            //rayon = K_FUNC::E_min(rayon, dx*dx+dy*dy+dz*dz);
          }
          rayon = rayon / nbV;
          //rayon = K_FUNC::E_max(1.e-10, rayon);
          
          // Essai de voir si le mouvement retourne la maille
          
          // Essaie de voir si le pt est sur le bord
          for (E_Int vi = 0; vi < nbV; vi++)
          {
            E_Int nov = v[vi]-1;
            dx = projx[nov]-cx[nov];
            dy = projy[nov]-cy[nov];
            dz = projz[nov]-cz[nov];
            distv = dx*dx + dy*dy + dz*dz;
            if (distv < delta*rayon) { edge = true; break; }
          }
          //printf("dist=%f rayon=%f edge=%d\n", dist, delta*rayon, edge);
        
          if (dist < delta*rayon && edge == true && movep[ind] >= 1.)
          {
            //alpha = pow(dist/(delta*rayon), 0.1);
            alpha = 0.0;
            f1x[ind] = cx[ind] + (1.-alpha)*(projx[ind]-cx[ind])+
              alpha*(f1x[ind]-cx[ind]);
            f1y[ind] = cy[ind] + (1.-alpha)*(projy[ind]-cy[ind])+
              alpha*(f1y[ind]-cy[ind]);
            f1z[ind] = cz[ind] + (1.-alpha)*(projz[ind]-cz[ind])+
              alpha*(f1z[ind]-cz[ind]);
            //printf("proj %f %f %f -> %f %f %f\n", cx[ind],cy[ind],cz[ind],f1x[ind],f1y[ind],f1z[ind]);
          }
          */
        }
      }
    }

    coord.setOneField(coordo, posx1, 1);
    coord.setOneField(coordo, posy1, 2);
    coord.setOneField(coordo, posz1, 3);
    nit++;
  }
}
