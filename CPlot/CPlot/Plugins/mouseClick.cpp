/*    
    Copyright 2013-2019 Onera.

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
#include "../cplot.h"
#include "../Data.h"

//=============================================================================
// Mouse plugins
//=============================================================================

int findNearestPoint(double xp, double yp, double zp,
                     StructZone* zone, int& ind, double& dist);
int findNearestPoint(double xp, double yp, double zp,
                     UnstructZone* zone, int& ind, double& dist);
int findElement(double xp, double yp, double zp,
                UnstructZone* zone, double& dist);
int findElement(double xp, double yp, double zp,
                StructZone* zone, double& dist);
int findFace(double xp, double yp, double zp, int elt, 
             UnstructZone* zone, double& dist);

//=============================================================================
// Click selects
//=============================================================================
void mouseClickSelect(Data* d, int button, int etat, int x, int y)
{
  d->dataMouseClickSelect(button, etat, x, y, 0, 0);
}
void mouseClickMultipleSelect(Data* d, int button, int etat, int x, int y)
{
  d->dataMouseClickSelect(button, etat, x, y, 1, 0);
}
void mouseClickAccurateSelect(Data* d, int button, int etat, int x, int y)
{
  d->dataMouseClickSelect(button, etat, x, y, 0, 1);
}
void mouseRightClickSelect(Data* d, int button, int etat, int x, int y)
{
  d->dataMouseRightClickSelect(button, etat, x, y);
}
void Data::dataMouseClickSelect(int button, int etat, int x, int y, 
                                int multiple, int accurate)
{
  if (etat == 1) {
    // Bouton relache
    return; }

  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLfloat winX, winY, winZ;
  GLdouble posX, posY, posZ;
  int n;

  //---------------------------------------------------------------
  // Compute mouse position: posX, posY, posZ in real coordinates
  //---------------------------------------------------------------
  // Ca marche super, a condition de cliquer sur un pixel!!
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);

  winX = (float)x;
  winY = (float)viewport[3] - (float)y;
  glReadPixels(x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
  // winZ: devrait etre corrige a cause du polygon offset + incertitude
  // du z-buffer
  if (winZ == 1.) return; // pas d'objet sous le click
  
  // Coord dans le monde reel
  gluUnProject(winX, winY, winZ, modelview, projection, viewport, 
               &posX, &posY, &posZ);

  ptrState->currentMousePosX = posX; ptrState->currentMousePosY = posY;
  ptrState->currentMousePosZ = posZ;

  // origine du rayon: on peut s'en servir pour calculer l'intersection
  // du rayon avec le maillage
  //double xorig, yorig, zorig;
  //gluUnProject( winX, winY, 0., modelview, projection, viewport, 
  //              &xorig, &yorig, &zorig);

  // Blocs contenant un point le plus proche de P
  int ret; int zone, ind, indE; double dist;
  ret = findBlockContaining(posX, posY, posZ, zone, ind, indE, dist);

  // Projection posX, posY, posZ sur la face adhoc
  //projectP(posX, posY, posZ, _zones[zone], indE);
  //printf("click: %f %f %f : %d\n", posX, posY, posZ, ret);

  // save selected status before click
  for (int i = 0; i < _numberOfZones; i++)
    _zones[i]->previouslySelected = _zones[i]->selected;

  Zone* z = _zones[zone];
  if (ret == 1 && (z->active == 1 || ptrState->ghostifyDeactivatedZones == 1))
  {
    if (accurate == 1)
    {
      ptrState->activePointX = z->x[ind];
      ptrState->activePointY = z->y[ind];
      ptrState->activePointZ = z->z[ind];
      ptrState->activePointZBuf = winZ;
    }
    else
    {
      // Il faudrait projeter P sur les facettes de l'element pour etre
      // plus precis
      ptrState->activePointX = posX;
      ptrState->activePointY = posY;
      ptrState->activePointZ = posZ;
      ptrState->activePointZBuf = winZ;
    }
    if (zone < _numberOfStructZones)
    {
      StructZone* zz = (StructZone*)z;
      int ni = zz->ni; 
      int nj = zz->nj;
      int k = ind / (ni*nj);
      int j = (ind - k*ni*nj)/ni;
      int i = ind - k*ni*nj - j*ni;
      ptrState->activePointI = i+1;
      ptrState->activePointJ = j+1;
      ptrState->activePointK = k+1;
    }
    else
    {
      ptrState->activePointI = ind; // indice du noeud le plus proche
      ptrState->activePointJ = indE; // indice de l'element contenant P
      UnstructZone* zz = (UnstructZone*)z;
      if (zz->eltType != 10) // autre que NGON
      {
        int* c = zz->connect;
        int size = zz->eltSize;
        int ne = zz->ne;
        int v = 0;
        for (v = 0; v < size; v++)
        {
          if (c[indE+v*ne] == ind+1) break;
        }
        ptrState->activePointK = -v-1;
      }
      else ptrState->activePointK = findFace(posX, posY, posZ, indE, zz, dist);
    }
    for (n = 0; n < z->nfield; n++)
    {
      double* f = z->f[n];
      ptrState->activePointF[n] = f[ind];
    }

    int previousSelected = ptrState->selectedZone;
    ptrState->selectedZone = zone+1;

    Zone* z = _zones[zone];

    if (ptrState->selectionStyle == 2) ptrState->ktimer = 5;

    if (multiple == 1)
    {
      if (z->selected == 1)
      {
        z->selected = 0;
        if (previousSelected == ptrState->selectedZone) ptrState->selectedZone = 0;
        else ptrState->selectedZone = previousSelected;
      }
      else z->selected = 1;
    }
    else // not multiple select
    {
      for (int i = 0; i < _numberOfZones; i++) _zones[i]->selected = 0;
      z->selected = 1;
    }
  }
}

//=============================================================================
// Shift + Right Click: deactivate zone
//=============================================================================
void Data::dataMouseRightClickSelect(int button, int etat, int x, int y)
{
  if (etat == 1) return;
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLfloat winX, winY, winZ;
  GLdouble posX, posY, posZ;

  glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
  glGetDoublev( GL_PROJECTION_MATRIX, projection );
  glGetIntegerv( GL_VIEWPORT, viewport );

  winX = (float)x;
  winY = (float)viewport[3] - (float)y;
  glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
  if (winZ == 1.) return; // pas d'objet sous le click

  gluUnProject( winX, winY, winZ, modelview, projection, viewport, 
                &posX, &posY, &posZ);
  // Des erreurs peuvent exister sur posX, posY, posZ.
  int ret; int zone, ind, indE; double dist;
  ret = findBlockContaining(posX, posY, posZ, zone, ind, indE, dist);
  
  Zone* z = _zones[zone];
  if (ret == 1 && (z->active == 1 || ptrState->ghostifyDeactivatedZones == 1))
  {
  
    if (z->active == 1)
    {
      z->active = 0;
      if (ptrState->deactivatedZones == NULL)
      {
        struct chain_int* ci;
        ci = (struct chain_int*)malloc(sizeof(struct chain_int));
        ci->value = zone+1;
        ci->next = NULL;
        ptrState->deactivatedZones = ci;
      }
      else
      {
        struct chain_int* ci = ptrState->deactivatedZones;
        while (ci->next != NULL) ci = ci->next;
        ci->next = (struct chain_int*)malloc(sizeof(struct chain_int));
        ci = ci->next;
        ci->value = zone+1;
        ci->next = NULL;
      }
    }
    else z->active = 1;
  }
}

//=============================================================================
// Retourne le no du bloc duquel P est le plus proche
// IN: x,y,z: coord. du point P a tester
// OUT: zone: numero du bloc le plus proche de P
// OUT: ind: indice du point le plus proche de P appartenant a zone
// OUT: indE: indice de l'element le plus proche de P appartenant a zone
// OUT: d: distance entre P et le point ind
// Retourne 1 si une zone a ete trouvee, 0 sinon.
//=============================================================================
int Data::findBlockContaining(double x, double y, double z,
                              int& zone, int& ind, int& indE,
                              double& dist)
{
  int nz, indl, inde;
  double xmi, ymi, zmi, xma, yma, zma;
  double d, dn, de;
  int nzmin = -1;
  double dmin = 1.e6; // distMin node/element
  double dminNode = 1.e6; // distMin node
  double dminElt = 1.e6; // distMin element
  ind = 0; indE = 0;
  double eps = 0.05;
  d = MAX(xmax-xmin, ymax-ymin);
  d = MAX(d, zmax-zmin);
  eps = d*eps;
  //printf("eps: %f\n", eps);

  for (nz = 0; nz < _numberOfZones; nz++)
  {
    Zone* zone = _zones[nz];
    if (zone->active == 1 || 
        (zone->active == 0 && ptrState->ghostifyDeactivatedZones == 1))
    {
      xma = zone->xmax + eps;
      xmi = zone->xmin - eps;
      yma = zone->ymax + eps;
      ymi = zone->ymin - eps;
      zma = zone->zmax + eps;
      zmi = zone->zmin - eps;
      //printf("try: zone %f %f %f %f %f %f\n", xmi, xma, ymi, yma, zmi, zma);
      if (x >= xmi && x <= xma && 
          y >= ymi && y <= yma && 
          z >= zmi && z <= zma)
      {        
        if (nz < _numberOfStructZones)
        {
          findNearestPoint(x, y, z, (StructZone*)zone, indl, dn);
          inde = findElement(x, y, z, (StructZone*)zone, de);
        }
        else 
        {
          findNearestPoint(x, y, z, (UnstructZone*)zone, indl, dn);
          inde = findElement(x, y, z, (UnstructZone*)zone, de);
        }
        d = MIN(dn, de);
        if (zone->active == 0) d = d + eps*0.01; // malus

        if (d < dmin-1.e-10)
        {
          dmin = d; nzmin = nz; ind = indl; indE = inde;
          dminElt = de; dminNode = dn;
        }
        else if (d <= dmin) // meme point mini: compare dn ou de
        {
          if (de < dminElt && _zones[nzmin]->dim != 1)
          {
            dmin = d; nzmin = nz; ind = indl; indE = inde; 
            dminElt = de; dminNode = dn;
          }
          else if (dn < dminNode && _zones[nzmin]->dim != 1)
          {
            dmin = d; nzmin = nz; ind = indl; indE = inde; 
            dminElt = de; dminNode = dn;
          }
        }
      }
    }
  }

  zone = nzmin;
  dist = dmin;
  if (nzmin == -1) return 0;
  else return 1;
}

//=============================================================================
// pour les grilles non structurees
//=============================================================================
int findNearestPoint(double xp, double yp, double zp,
                     UnstructZone* zone, int& ind, double& dist)
{
  double d, dx, dy, dz;
  dist = 1.e6;
  int npts = zone->npts;
  double* x = zone->x;
  double* y = zone->y;
  double* z = zone->z;
  ind = 0;
 
  for (int i = 0; i < npts; i++)
  {
    dx = x[i]-xp; dy = y[i]-yp; dz = z[i]-zp;
    d = dx*dx + dy*dy + dz*dz;
    if (d < dist)
    {
      dist = d; ind = i;
    }
  }
  return 0;
}
//=============================================================================
// pour les grilles structurees
//=============================================================================
int findNearestPoint(double xp, double yp, double zp,
                     StructZone* zone, int& ind, double& dist)
{
  double d, dx, dy, dz;
  dist = 1.e6;
  double* x = zone->x;
  double* y = zone->y;
  double* z = zone->z;
  int iplane = zone->iPlane; // from 0 to ni-1, -1 means min-max planes
  int jplane = zone->jPlane;
  int kplane = zone->kPlane;
  int ni = zone->ni;
  int nj = zone->nj;
  int nk = zone->nk;
  int ninj = ni*nj;
  ind = 0;
  int inds, inc; 

  if (iplane >= 0) // un plan
  {
    inc = iplane;
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
      {
        inds = inc + j*ni + k* ninj;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
      }
  }
  else // plans imin,imax
  {
    inc = ni-1;
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
      {
        inds = j*ni + k* ninj;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
        inds = inds + inc;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
      }
  }
  if (jplane >= 0) 
  {
    inc = jplane*ni;
    for (int k = 0; k < nk; k++)
      for (int i = 0; i < ni; i++)
      {
        inds = i + inc + k* ninj;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
      }
  }
  else 
  {
    inc = (nj-1)*ni;
    for (int k = 0; k < nk; k++)
      for (int i = 0; i < ni; i++)
      {
        inds = i +  k*ninj;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
        
        inds = inds + inc;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
      }
  }
  if (kplane >= 0) 
  {
    inc = kplane*ninj;
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
      {
        inds = i + j*ni + inc;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
      }
  } 
  else // plans kmin/kmax
  {
    inc = (nk-1)*ninj;
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
      {
        inds = i + j*ni;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
        inds = inds + inc;
        dx = x[inds]-xp; dy = y[inds]-yp; dz = z[inds]-zp;
        d = dx*dx + dy*dy + dz*dz;
        if (d < dist)
        { dist = d; ind = inds;}
      }
  }
  return 0;
}
//=============================================================================
// Trouve l'element dont le centre est le plus proche de P (xp, yp, zp). 
// (pour une zone non structuree)
// Retourne distMin: la distance du centre a P
//=============================================================================
int findElement(double xp, double yp, double zp,
                UnstructZone* zone, double& distMin)
{
  int ne = zone->ne;
  int nv = zone->eltSize;
  double nvi = 1./MAX(nv,1);
  int* c = zone->connect;
  double* x = zone->x;
  double* y = zone->y;
  double* z = zone->z;
  int i, indl, v;
  double xc, yc, zc, dist, dx, dy, dz;
  distMin = 1.e6; 
  int best = 0;
  if (zone->eltType != 10) // basic elements
  {
    for (i = 0; i < ne; i++)
    {
      xc = 0.; yc = 0.; zc = 0.;
      for (v = 0; v < nv; v++)
      {
        indl = c[i+ne*v]-1;
        xc += x[indl]; yc += y[indl]; zc += z[indl];
      }
      xc = xc * nvi; yc = yc * nvi; zc = zc * nvi;
      dx = xp-xc; dy = yp-yc; dz = zp-zc;
      dist = dx*dx + dy*dy + dz*dz;
      if (dist < distMin)
      {distMin = dist; best = i;}
    }
  }
  else
  { // NGONS
    int* ptr = PTRELTS(c); // ptr sur les elts
    int nf, f, p, np, rt;
    int* lp;
    int* posf = zone->posFaces;
    for (i = 0; i < ne; i++)
    {
      xc = 0.; yc = 0.; zc = 0.; rt = 0;
      nf = ptr[0];
      for (f = 0; f < nf; f++) // pour chaque face
      {
        lp = &c[posf[ptr[f+1]-1]];
        np = lp[0];
        for (p = 0; p < np; p++)
        {
          xc += x[lp[p+1]-1];
          yc += y[lp[p+1]-1];
          zc += z[lp[p+1]-1]; rt++;
        }
      }
      ptr += nf+1;
      rt = MAX(rt, 1);
      xc = xc * 1./rt; yc = yc * 1./rt; zc = zc * 1./rt;
      dx = xp-xc; dy = yp-yc; dz = zp-zc;
      dist = dx*dx + dy*dy + dz*dz;
      if (dist < distMin)
      {distMin = dist; best = i;}
    }
  }
  return best;
}

//=============================================================================
// Trouve l'element dont le centre est le plus proche de P (xp, yp, zp). 
// (pour une zone structuree)
// Retourne: distMin: la distance de P au centre.
//=============================================================================
int findElement(double xp, double yp, double zp,
                StructZone* zone, double& distMin)
{
  int ni = zone->ni;
  int nj = zone->nj;
  int nk = zone->nk;
  int ninj = ni*nj;
  double* x = zone->x;
  double* y = zone->y;
  double* z = zone->z;
  int iplane = zone->iPlane; // from 0 to ni-1, -1 means min-max planes
  int jplane = zone->jPlane;
  int kplane = zone->kPlane;
  double xc, yc, zc, dist;
  int indElt = 0;
  int inc, incp, indE; 
  int ii1, jj1, kk1, ii2, jj2, kk2, ind1, ind2, ind3, ind4, ind5;
  int ind6, ind7, ind8;
  int isup, jsup, ksup;
  int ni1 = ni-1, nj1 = nj-1, nk1 = nk-1, ni1nj1 = ni1*nj1;
  isup = MAX(ni1,1); jsup = MAX(nj1,1); ksup = MAX(nk1,1);
  distMin = 1.e6;

  if (iplane >= 0) // un plan
  {
    ii1 = iplane; ii2 = MIN(ii1+1, ni1);
    for (int ks = 0; ks < ksup; ks++)
    {
      kk1 = ks; kk2 = MIN(kk1+1, nk1);
      for (int js = 0; js < jsup; js++)
      {
        jj1 = js; jj2 = MIN(jj1+1, nj1);
        indE = ii1 + jj1*ni1 + kk1*ni1nj1;
        ind1 = ii1 + jj1*ni + kk1*ninj;
        ind2 = ii2 + jj1*ni + kk1*ninj;
        ind3 = ii1 + jj2*ni + kk1*ninj;
        ind4 = ii2 + jj2*ni + kk1*ninj;
        ind5 = ii1 + jj1*ni + kk2*ninj;
        ind6 = ii2 + jj1*ni + kk2*ninj;
        ind7 = ii1 + jj2*ni + kk2*ninj;
        ind8 = ii2 + jj2*ni + kk2*ninj;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
      } 
    }
  }
  else // plans imin,imax
  {
    inc = MAX(ni1-1,0);
    if (ni == 1) incp = 0;
    else incp = 1;
    for (int ks = 0; ks < ksup; ks++)
    {
      kk1 = ks; kk2 = MIN(kk1+1, nk1);
      for (int js = 0; js < jsup; js++)
      {
        jj1 = js; jj2 = MIN(jj1+1, nj1);
        indE = jj1*ni1 + kk1*ni1nj1;
        ind1 = jj1*ni + kk1*ninj; ind2 = ind1+incp;
        ind3 = jj2*ni + kk1*ninj; ind4 = ind3+incp;
        ind5 = jj1*ni + kk2*ninj; ind6 = ind5+incp;
        ind7 = jj2*ni + kk2*ninj; ind8 = ind7+incp;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
        // plan imax
        indE = indE + inc;
        ind1 = ind1+inc; ind2 = ind2+inc; ind3 = ind3+inc; ind4 = ind4+inc;
        ind5 = ind5+inc; ind6 = ind6+inc; ind7 = ind7+inc; ind8 = ind8+inc;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
      }
    }
  }
  if (jplane >= 0) // un plan
  {
    jj1 = jplane; jj2 = MIN(jj1+1, nj1);
    for (int ks = 0; ks < ksup; ks++)
    {
      kk1 = ks; kk2 = MIN(kk1+1, nk1); 
      for (int is = 0; is < isup; is++)
      {
        ii1 = is; ii2 = MIN(ii1+1, ni1);
        indE = ii1 + jj1*ni1 + kk1*ni1nj1;
        ind1 = ii1 + jj1*ni + kk1*ninj;
        ind2 = ii2 + jj1*ni + kk1*ninj;
        ind3 = ii1 + jj2*ni + kk1*ninj;
        ind4 = ii2 + jj2*ni + kk1*ninj;
        ind5 = ii1 + jj1*ni + kk2*ninj;
        ind6 = ii2 + jj1*ni + kk2*ninj;
        ind7 = ii1 + jj2*ni + kk2*ninj;
        ind8 = ii2 + jj2*ni + kk2*ninj;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
      }
    }
  }
  else // plans jmin et jmax
  {
    inc = MAX((nj1-1)*ni,0);
    if (nj == 1) incp = 0;
    else incp = ni;
    for (int ks = 0; ks < ksup; ks++)
    {
      kk1 = ks; kk2 = MIN(kk1+1, nk1);
      for (int is = 0; is < isup; is++)
      {
        ii1 = is; ii2 = MIN(ii1+1, ni1);
        indE = ii1 + kk1*ni1nj1;
        ind1 = ii1 + kk1*ninj; ind3 = ind1+incp;
        ind2 = ii2 + kk1*ninj; ind4 = ind2+incp;
        ind5 = ii1 + kk2*ninj; ind7 = ind5+incp; 
        ind6 = ii2 + kk2*ninj; ind8 = ind6+incp;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
        // plan jmax
        indE = indE + MAX((nj1-1)*ni1,0);
        ind1 = ind1+inc; ind2 = ind2+inc; ind3 = ind3+inc; ind4 = ind4+inc;
        ind5 = ind5+inc; ind6 = ind6+inc; ind7 = ind7+inc; ind8 = ind8+inc;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
      }
    }
  }
  if (kplane >= 0) // un plan
  {
    kk1 = kplane; kk2 = MIN(kk1+1, nk1);
    for (int js = 0; js < jsup; js++)
    {
      jj1 = js; jj2 = MIN(jj1+1, nj1);
      for (int is = 0; is < isup; is++)
      {
        ii1 = is; ii2 = MIN(ii1+1, ni1);
        indE = ii1 + jj1*ni1 + kk1*ni1nj1;
        ind1 = ii1 + jj1*ni + kk1*ninj;
        ind2 = ii2 + jj1*ni + kk1*ninj;
        ind3 = ii1 + jj2*ni + kk1*ninj;
        ind4 = ii2 + jj2*ni + kk1*ninj;
        ind5 = ii1 + jj1*ni + kk2*ninj;
        ind6 = ii2 + jj1*ni + kk2*ninj;
        ind7 = ii1 + jj2*ni + kk2*ninj;
        ind8 = ii2 + jj2*ni + kk2*ninj;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
      }
    }
  }
  else // plans kmin et kmax
  {
    inc = MAX((nk1-1)*ninj, 0);
    if (nk == 1) incp = 0;
    else incp = ninj;
    for (int js = 0; js < jsup; js++)
    {
      jj1 = js; jj2 = MIN(jj1+1, nj1);
      for (int is = 0; is < isup; is++)
      {
        ii1 = is; ii2 = MIN(ii1+1, ni1);
        indE = ii1 + jj1*ni1;
        ind1 = ii1 + jj1*ni; ind5 = ind1 + incp;
        ind2 = ii2 + jj1*ni; ind6 = ind2 + incp;
        ind3 = ii1 + jj2*ni; ind7 = ind3 + incp;
        ind4 = ii2 + jj2*ni; ind8 = ind4 + incp;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
        // plan kmax
        indE = indE + MAX((nk1-1)*ni1nj1,0);
        ind1 = ind1+inc; ind2 = ind2+inc; ind3 = ind3+inc; ind4 = ind4+inc;
        ind5 = ind5+inc; ind6 = ind6+inc; ind7 = ind7+inc; ind8 = ind8+inc;
        xc = x[ind1]+x[ind2]+x[ind3]+x[ind4]+x[ind5]+x[ind6]+x[ind7]+x[ind8];
        yc = y[ind1]+y[ind2]+y[ind3]+y[ind4]+y[ind5]+y[ind6]+y[ind7]+y[ind8];
        zc = z[ind1]+z[ind2]+z[ind3]+z[ind4]+z[ind5]+z[ind6]+z[ind7]+z[ind8];
        xc = 0.125*xc; yc = 0.125*yc; zc = 0.125*zc;
        dist = (xp-xc)*(xp-xc) + (yp-yc)*(yp-yc) + (zp-zc)*(zp-zc);
        if (dist < distMin)
        {indElt = indE; distMin = dist;}
      }
    }
  }
  return indElt;
}

//=============================================================================
// NGON: Trouve la face de l'element elt dont le centre est le plus
// proche de xp,yp,zp
// le no de la face commence a 0, elt commence a 0
//=============================================================================
int findFace(double xp, double yp, double zp, int elt, 
             UnstructZone* zone, double& dist)
{
  int* c = zone->connect;
  int* ptr = PTRELTS(c); // ptr sur les elts
  int f, p, np, rt, nf, i, ne;
  ne = zone->ne;
  double xc, yc, zc, d;
  double* x = zone->x;
  double* y = zone->y;
  double* z = zone->z;
  int* lp;
  int* posf = zone->posFaces;
  dist = 1.e6;
  int face = 0;

  for (i = 0; i < ne; i++)
  {
    nf = ptr[0];
    if (i == elt)
    {
      for (f = 0; f < nf; f++) // pour chaque face
      {
        xc = 0.; yc = 0.; zc = 0.;
        lp = &c[posf[ptr[f+1]-1]];
        np = lp[0];
        for (p = 0; p < np; p++)
        {
          xc += x[lp[p+1]-1];
          yc += y[lp[p+1]-1];
          zc += z[lp[p+1]-1];
        }
        rt = MAX(np, 1);
        xc = xc / rt; yc = yc / rt; zc = zc / rt;
        d = (xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp);
        if (d < dist)
        { dist = d; face = ptr[f+1]-1; }
      }
      break;
    }
    ptr += nf+1;
  }
  return face;
}
