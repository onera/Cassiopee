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
#include "../Data.h"

//=============================================================================
/* 
   Display a structured mesh edges.
   Affiche les lignes imin, imax, jmin, jmax, kmin, kmax.
*/
//=============================================================================
void Data::displaySEdges()
{
  if (_numberOfStructZones == 0) return;
  if (ptrState->edgifyDeactivatedZones == 0 && ptrState->edgifyActivatedZones == 0)
    return;
  int zone;

  glColor4f(0.8, 0.8, 0.8, 1.);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  zone = 0;
  while (zone < _numberOfStructZones)
  {
    StructZone* zonep = _szones[zone];
    if (isInFrustum(zonep, _view) == 1)
    {
      // if zone is active and in frustum
      if ((zonep->active == 1 && ptrState->edgifyActivatedZones == 1) ||
          (zonep->active == 0 && ptrState->edgifyDeactivatedZones == 1))
        
      {
        int ni, nj, nk, n1;
        ni = zonep->ni; nj = zonep->nj; nk = zonep->nk;
        if (ptrState->dim == 2) nk = 1;
        int ni1 = ni-1; int nj1 = nj-1; int nk1 = nk-1; 
        int nij = ni*nj;
        double* x = zonep->x; double* y = zonep->y; double* z = zonep->z;

        // ligne 1
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < ni; i++) 
        { n1 = i; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 2
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < ni; i++) 
        { n1 = i+nj1*ni; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 3
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < ni; i++) 
        { n1 = i+nk1*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 4
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < ni; i++) 
        { n1 = i+nj1*ni+nk1*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 5
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < nj; j++) 
        { n1 = j*ni; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();
        
        // ligne 6
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < nj; j++) 
        { n1 = ni1+j*ni; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();
        
        // ligne 7
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < nj; j++) 
        { n1 = j*ni+nk1*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 8
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < nj; j++) 
        { n1 = ni1+j*ni+nk1*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 9
        glBegin(GL_LINE_STRIP);
        for (int k = 0; k < nk; k++) 
        { n1 = k*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 10
        glBegin(GL_LINE_STRIP);
        for (int k = 0; k < nk; k++) 
        { n1 = ni1+k*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 11
        glBegin(GL_LINE_STRIP);
        for (int k = 0; k < nk; k++) 
        { n1 = nj1*ni+k*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();

        // ligne 12
        glBegin(GL_LINE_STRIP);
        for (int k = 0; k < nk; k++) 
        { n1 = ni1+nj1*ni+k*nij; glVertex3d(x[n1], y[n1], z[n1]); }
        glEnd();
      }
    }
    zone++;
  }
 
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_BLEND);
  ptrState->alpha = 1.; 
  glColor3f(1., 1., 1.);
}

