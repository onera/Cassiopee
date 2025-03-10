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
#include "../Data.h"

//============================================================================
/* 
   Lighten the scene 
   type=0: hard, one-sided
   type=1: soft, one-sided
   type=2: hard, two-sided
   type=3: soft, two-sided
*/
//============================================================================
void Data::light(E_Int type)
{
  GLfloat ambient[4];
  GLfloat diffuse[4];
  GLfloat specular[4];
  GLfloat shininess[1];
  switch (type)
  {
    case 0: // hard
    case 2:
      ambient[0] = 0.04; // old: 0.01
      ambient[1] = 0.04;
      ambient[2] = 0.04; 
      ambient[3] = 1.0;
      diffuse[0] = 1.;
      diffuse[1] = 1.;
      diffuse[2] = 1.;
      diffuse[3] = 1.;
      specular[0] = 0.3; // old: 1
      specular[1] = 0.3;
      specular[2] = 0.3;
      specular[3] = 0.3;
      shininess[0] = 50; // old: 60
      break;

    case 1: // soft
    case 3:
      ambient[0] = 0.1; // old: 0.2
      ambient[1] = 0.1;
      ambient[2] = 0.1; 
      ambient[3] = 1.;
      diffuse[0] = 1.;
      diffuse[1] = 1.;
      diffuse[2] = 1.;
      diffuse[3] = 1.;
      specular[0] = 0.2;
      specular[1] = 0.2;
      specular[2] = 0.2;
      specular[3] = 0.2;
      shininess[0] = 80;
      break;
  }
  
  // Enable calculation of lighting
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  // 0.0 indique une lumiere directionelle
  // x,y,z sont la direction
#ifdef __SHADERS__
  // Position fixe dans le repere du modele (pour le shader)
  //GLfloat position[] = { 0., 0., 0., 0.0 };
  double d1x = _view.xcam-_view.xeye;
  double d1y = _view.ycam-_view.yeye;
  double d1z = _view.zcam-_view.zeye;
  double d2x = _view.dirx;
  double d2y = _view.diry;
  double d2z = _view.dirz;
  double d3x = d1y*d2z-d1z*d2y;
  double d3y = d1z*d2x-d1x*d2z;
  double d3z = d1x*d2y-d1y*d2x;
  double r = d1x*d1x+d1y*d1y+d1z*d1z;
  double n = d3x*d3x+d3y*d3y+d3z*d3z;
  n = 1./fmax(sqrt(n), 1.e-12);
  d3x = d3x*n; d3y=d3y*n; d3z=d3z*n;
  n = d2x*d2x+d2y*d2y+d2z*d2z;
  n = 1./fmax(sqrt(n), 1.e-12);
  d2x = d2x*n; d2y=d2y*n; d2z=d2z*n;
  r = sqrt(r)*3.14;
  float xl = ptrState->lightOffsetX * r * d3x + ptrState->lightOffsetY * r * d2x;
  float yl = ptrState->lightOffsetX * r * d3y + ptrState->lightOffsetY * r * d2y;
  float zl = ptrState->lightOffsetX * r * d3z + ptrState->lightOffsetY * r * d2z;

  GLfloat position[] = { xl, yl, zl, 0.0 };
#else
  // Position dependant de la BB
  //GLfloat position[] = { (xmin+xmax)*0.5, ymin-2., zmax+2., 0.0 };
  // Eclairage a partir de l'observateur
  GLfloat position[] = { _view.xcam, _view.ycam, _view.zcam, 0.0 };
#endif

  // Define material of objects
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_BACK, GL_SPECULAR, specular);
  glMaterialfv(GL_BACK, GL_SHININESS, shininess);
  glMaterialfv(GL_BACK, GL_DIFFUSE, diffuse);

  // Define colors and position of LIGHT0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT0, GL_POSITION, position);

  // Enable different colors material for light
  //glEnable(GL_COLOR_MATERIAL);
  //glColorMaterial(GL_FRONT, GL_DIFFUSE);
  //glColorMaterial(GL_BACK, GL_DIFFUSE);

  if (type == 2 || type == 3)
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
              
}

//=============================================================================
void Data::noLight()
{
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  glDisable(GL_COLOR_MATERIAL);
}
