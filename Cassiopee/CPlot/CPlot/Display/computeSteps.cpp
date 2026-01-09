/*    
    Copyright 2013-2026 ONERA.

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
// Calcul les steps en se basant sur la distance a la BBOX
//=============================================================================
void Data::computeSteps0(StructZone* zonep, 
                         E_Int& stepi, E_Int& stepj, E_Int& stepk)
{
  double Di = zonep->Di;
  double Dj = zonep->Dj;
  double Dk = zonep->Dk;

  // Step setting
  double d = dist2BB(_view.xcam, _view.ycam, _view.zcam,
                     zonep->xmin, zonep->ymin, zonep->zmin,
                     zonep->xmax, zonep->ymax, zonep->zmax);
  
  double alpha = (1 + _pref.speed)*1.e-4;
  double beta = 0.5;
  d = pow(d, beta);
  d = d*alpha;

  stepi = (E_Int)(d*Di);
  stepj = (E_Int)(d*Dj);
  stepk = (E_Int)(d*Dk);
  stepi = MAX(E_Int(1), stepi);
  stepj = MAX(E_Int(1), stepj);
  stepk = MAX(E_Int(1), stepk);
  stepi = MIN(stepi, E_Int(20));
  stepj = MIN(stepj, E_Int(20));
}

//=============================================================================
// Calcul les steps en ses basant sur la distance a la BBOX
//=============================================================================
void Data::computeSteps1(StructZone* zonep, 
                         E_Int& stepi, E_Int& stepj, E_Int& stepk)
{
  double Di = zonep->Di;
  double Dj = zonep->Dj;
  double Dk = zonep->Dk;

  // Step setting
  double d = dist2BB(_view.xcam, _view.ycam, _view.zcam,
                     zonep->xmin, zonep->ymin, zonep->zmin,
                     zonep->xmax, zonep->ymax, zonep->zmax);

  // Code comme les mouvements de souris
  d = (sqrt(d)/dmoy)*2.;
  //printf("dist %f %f\n", d, Di);

  stepi = (E_Int)(d*Di);
  stepj = (E_Int)(d*Dj);
  stepk = (E_Int)(d*Dk);
  stepi = MAX(E_Int(1), stepi);
  stepj = MAX(E_Int(1), stepj);
  stepk = MAX(E_Int(1), stepk);
  stepi = MIN(stepi, E_Int(20));
  stepj = MIN(stepj, E_Int(20));
  stepk = MIN(stepk, E_Int(20));
}

//=============================================================================
// Calcul les steps en se basant sur les pixels
// On calcule les coordonnees d'un edge moyen dans les coord ecran
// Si il fait moins d'un pixel, on augmente le step
//=============================================================================
void Data::computeSteps(StructZone* zonep, 
                        E_Int& stepi, E_Int& stepj, E_Int& stepk)
{
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLdouble winX, winY, winZ;
  GLdouble posX, posY, posZ;
  GLdouble winX2, winY2, winZ2;
  GLdouble posX2, posY2, posZ2;
  GLdouble winX3, winY3, winZ3;
  GLdouble posX3, posY3, posZ3;
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);
  posX = zonep->x[0];
  posY = zonep->y[0];
  posZ = zonep->z[0];
  gluProject(posX, posY, posZ, modelview, projection, viewport, 
             &winX, &winY, &winZ);
  double Di = 1./zonep->Di;
  posX2 = zonep->x[0] + Di;
  posY2 = zonep->y[0] + Di;
  posZ2 = zonep->z[0] + Di;
  gluProject(posX2, posY2, posZ2, modelview, projection, viewport, 
             &winX2, &winY2, &winZ2);
  double Dx = K_FUNC::E_abs(winX-winX2);
  double Dy = K_FUNC::E_abs(winY-winY2);
  double Dz = K_FUNC::E_abs(winZ-winZ2);
  double DD;
  DD = MAX(Dx, Dy); DD = MAX(DD, Dz);
  posX3 = zonep->x[0] + 0.;
  posY3 = zonep->y[0] + Di;
  posZ3 = zonep->z[0] + 0;
  gluProject(posX3, posY3, posZ3, modelview, projection, viewport, 
             &winX3, &winY3, &winZ3);
  Dx = K_FUNC::E_abs(winX-winX3);
  Dy = K_FUNC::E_abs(winY-winY3);
  Dz = K_FUNC::E_abs(winZ-winZ3);
  DD = MAX(DD, Dx); DD = MAX(DD, Dy); DD = MAX(DD, Dz);
  //printf("DD=%f \n", DD);
  stepi = E_Int(0.8/DD)+1; 
  stepi = MIN(stepi, E_Int(20));
  stepi = MAX(stepi, 1);
  stepj = stepi; stepk = stepi;
  //if (stepi >= 2) printf("%d\n", stepi);
}
