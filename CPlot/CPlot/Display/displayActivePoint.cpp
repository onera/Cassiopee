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

//============================================================================
/*
  Display active point 
  - red box.
  - pour les zones non-structurees: axes x,y,z
  - pour les zones Structurees: directions i,j,k
  - optionel: HUD
*/
//============================================================================
void Data::displayActivePoint()
{
  int nz = ptrState->selectedZone;
  if (nz == 0) return;
  if (_zones[nz-1]->active == 0) return;

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(-1., -10.);
  if (_view.clipping >= 1) glDisable(GL_DEPTH_TEST);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1., 0.9, 0.9, 0.4);
  double x, y, z, dx, dy, dz, r, p;
  double ex, ey, ez;

  x = ptrState->activePointX;
  y = ptrState->activePointY;
  z = ptrState->activePointZ;

  ex = _view.xcam - _view.xeye;
  ey = _view.ycam - _view.yeye;
  ez = _view.zcam - _view.zeye;  

  p = ex*ex+ey*ey+ez*ez;
  p = 1./K_FUNC::E_max(sqrt(p), 1.e-12); // norme de e

  dx = _view.xcam - x;
  dy = _view.ycam - y;
  dz = _view.zcam - z;

  r = ex*dx + ey*dy + ez*dz;
  dx = r*ex*p; dy = r*ey*p; dz = r*ez*p;

  //r = dx*dx+dy*dy+dz*dz;
  r = sqrt(r)*0.005;
  p = 3.5*r;

  /* box */
  glBegin(GL_QUADS);
  glNormal3f(0,0,-1);
  glVertex3d(x-r, y-r, z-r);
  glVertex3d(x-r, y+r, z-r);
  glVertex3d(x+r, y+r, z-r);
  glVertex3d(x+r, y-r, z-r);
  
  glNormal3f(0,0,1);
  glVertex3d(x-r, y-r, z+r);
  glVertex3d(x-r, y+r, z+r);
  glVertex3d(x+r, y+r, z+r);
  glVertex3d(x+r, y-r, z+r);
  
  glNormal3f(-1,0,0);
  glVertex3d(x-r, y-r, z-r);
  glVertex3d(x-r, y+r, z-r);
  glVertex3d(x-r, y+r, z+r);
  glVertex3d(x-r, y-r, z+r);

  glNormal3f(1,0,0);
  glVertex3d(x+r, y-r, z-r);
  glVertex3d(x+r, y+r, z-r);
  glVertex3d(x+r, y+r, z+r);
  glVertex3d(x+r, y-r, z+r);

  glNormal3f(0,-1,0);
  glVertex3d(x-r, y-r, z-r);
  glVertex3d(x+r, y-r, z-r);
  glVertex3d(x+r, y-r, z+r);
  glVertex3d(x-r, y-r, z+r);

  glNormal3f(0,1,0);
  glVertex3d(x-r, y+r, z-r);
  glVertex3d(x+r, y+r, z-r);
  glVertex3d(x+r, y+r, z+r);
  glVertex3d(x-r, y+r, z+r);

  glEnd();

  nz = nz-1;
  
  if (nz < _numberOfStructZones)
  {
    glColor4f(1., 0., 0.,1.);
    glLineWidth(4.);
    glBegin(GL_LINES);

    /* Structurees: axes */
    StructZone* zp = (StructZone*)_zones[nz];
    int ni = zp->ni;
    int nj = zp->nj;
    int nk = zp->nk;
    int nij = ni*nj;
    int is = MAX(ptrState->activePointI-1,0); is = MIN(is, ni-1);
    int js = MAX(ptrState->activePointJ-1,0); js = MIN(js, nj-1);
    int ks = MAX(ptrState->activePointK-1,0); ks = MIN(ks, nk-1);
    int ind = is + js*ni + ks*nij;
    double* px = zp->x;
    double* py = zp->y;
    double* pz = zp->z;
    double vx, vy, vz, n;
    int indp, ip, jp, kp;

    // axe i
    if (is < ni-1 && ni > 1)
    {
      ip = is+1; jp = js; kp = ks;
      indp = ip + jp*ni + kp*nij;
      vx = px[indp]-px[ind];
      vy = py[indp]-py[ind];
      vz = pz[indp]-pz[ind];
      n = sqrt(vx*vx + vy*vy + vz*vz);
      n = 1./K_FUNC::E_max(n, 1.e-10);
      vx *= n; vy *= n; vz *= n;
      vx = vx*p; vy = vy*p; vz = vz*p;
      glColor4f(1.,0.,0.,1.);
      glVertex3d(x, y, z);
      glVertex3d(x+vx, y+vy, z+vz);
    }
    else if (ni > 1)
    {
      ip = is-1; jp = js; kp = ks;
      indp = ip + jp*ni + kp*nij;
      vx = px[indp]-px[ind];
      vy = py[indp]-py[ind];
      vz = pz[indp]-pz[ind];
      n = sqrt(vx*vx + vy*vy + vz*vz);
      n = 1./K_FUNC::E_max(n, 1.e-10);
      vx *= n; vy *= n; vz *= n;
      vx = -vx*p; vy = -vy*p; vz = -vz*p;
      glColor4f(1.,0.,0.,1.);
      glVertex3d(x, y, z);
      glVertex3d(x+vx, y+vy, z+vz);
    }
    
    // axe j
    if (js < nj-1 && nj > 1)
    {
      ip = is; jp = js+1; kp = ks;
      indp = ip + jp*ni + kp*nij;
      vx = px[indp]-px[ind];
      vy = py[indp]-py[ind];
      vz = pz[indp]-pz[ind];
      n = sqrt(vx*vx + vy*vy + vz*vz);
      n = 1./K_FUNC::E_max(n, 1.e-10);
      vx *= n; vy *= n; vz *= n;
      vx = vx*p; vy = vy*p; vz = vz*p;
      glColor4f(0.,1.,0.,1.);
      glVertex3d(x, y, z);
      glVertex3d(x+vx, y+vy, z+vz);
    }
    else if (nj > 1)
    {
      ip = is; jp = js-1; kp = ks;
      indp = ip + jp*ni + kp*nij;
      vx = px[indp]-px[ind];
      vy = py[indp]-py[ind];
      vz = pz[indp]-pz[ind];
      n = sqrt(vx*vx + vy*vy + vz*vz);
      n = 1./K_FUNC::E_max(n, 1.e-10);
      vx *= n; vy *= n; vz *= n;
      vx = -vx*p; vy = -vy*p; vz = -vz*p;
      glColor4f(0.,1.,0.,1.);
      glVertex3d(x, y, z);
      glVertex3d(x+vx, y+vy, z+vz);
    }

    // axe k
    if (ks < nk-1 && nk > 1)
    {
      ip = is; jp = js; kp = ks+1;
      indp = ip + jp*ni + kp*nij;
      vx = px[indp]-px[ind];
      vy = py[indp]-py[ind];
      vz = pz[indp]-pz[ind];
      n = sqrt(vx*vx + vy*vy + vz*vz);
      n = 1./K_FUNC::E_max(n, 1.e-10);
      vx *= n; vy *= n; vz *= n;
      vx = vx*p; vy = vy*p; vz = vz*p;
      glColor4f(0.,0.,1.,1.);
      glVertex3d(x, y, z);
      glVertex3d(x+vx, y+vy, z+vz);
    }
    else if (nk > 1)
    {
      ip = is; jp = js; kp = ks-1;
      indp = ip + jp*ni + kp*nij;
      vx = px[indp]-px[ind];
      vy = py[indp]-py[ind];
      vz = pz[indp]-pz[ind];
      n = sqrt(vx*vx + vy*vy + vz*vz);
      n = 1./K_FUNC::E_max(n, 1.e-10);
      vx *= n; vy *= n; vz *= n;
      vx = -vx*p; vy = -vy*p; vz = -vz*p;
      glColor4f(1.,0.,1.,1.);
      glVertex3d(x, y, z);
      glVertex3d(x+vx, y+vy, z+vz);
    }
    glEnd();
  }
  else
  {
    /* Non structurees: axes */
    glLineWidth(4.);
    glColor4f(1., 0., 0.,1.);
    glBegin(GL_LINES);
    glVertex3d(x-p, y, z);
    glVertex3d(x+p, y, z);
    glVertex3d(x, y-p, z);
    glVertex3d(x, y+p, z);
    glVertex3d(x, y, z-p);
    glVertex3d(x, y, z+p);
    glEnd();
  }

  /* HUD (try) */
  glDisable(GL_DEPTH_TEST);
  double dirx, diry, dirz, vx, vy, vz, n;
  dirx = _view.dirx; diry = _view.diry; dirz = _view.dirz;
  vx = diry*dz - dirz*dy;
  vy = dirz*dx - dirx*dz;
  vz = dirx*dy - diry*dx;
  n = sqrt(vx*vx+vy*vy+vz*vz);
  n = 1./K_FUNC::E_max(n, 1.e-12);
  vx *= n; vy *= n; vz *= n;

  glLineWidth(2.);
  glColor4f(1., 0., 0.,1.);
  glBegin(GL_LINES);
  glVertex3d(x, y, z);
  glVertex3d(x+5*r*(vx+dirx), y+5*r*(vy+diry), z+5*r*(vz+dirz));
  glVertex3d(x+5*r*(vx+dirx), y+5*r*(vy+diry), z+5*r*(vz+dirz));
  glVertex3d(x+5*r*(2*vx+dirx), y+5*r*(2*vy+diry), z+5*r*(2*vz+dirz));
  glEnd();
  
  char msg[512]; char msg2[512]; msg2[0] = '\0';
  if (ptrState->selectedZone > 0)
  {
    Zone* zone = _zones[ptrState->selectedZone-1];
    if (ptrState->selectedZone-1 < _numberOfStructZones)
    {
      /*
      if (ptrState->mode == MESH) // xyz + ijk
      {
        sprintf(msg, "%g, %g, %g", 
                ptrState->activePointX, ptrState->activePointY, ptrState->activePointZ);
        sprintf(msg2, "(i=%d, j=%d, k=%d)",
                ptrState->activePointI, ptrState->activePointJ, abs(ptrState->activePointK));
      }*/
      if (ptrState->mode == SCALARFIELD)
      {
        sprintf(msg, "%s=%g",
                zone->varnames[ptrState->scalarField], 
                ptrState->activePointF[ptrState->scalarField]);
      }
      else if (ptrState->mode == VECTORFIELD)
      {
        sprintf(msg, "%s=%g, %s=%g, %s=%g",
                zone->varnames[ptrState->vectorField1], 
                zone->f[ptrState->vectorField1][ptrState->activePointI],
                zone->varnames[ptrState->vectorField2],
                zone->f[ptrState->vectorField2][ptrState->activePointI],
                zone->varnames[ptrState->vectorField3],
                zone->f[ptrState->vectorField3][ptrState->activePointI]);
      }
      else
        sprintf(msg,"%s (STRUCT)", zone->zoneName);
    }
    else
    {
      UnstructZone* zu = (UnstructZone*)zone;
      /*
      if (ptrState->mode == MESH) // xyz + np,ne
      {
        sprintf(msg, "%g, %g, %g", 
        ptrState->activePointX, ptrState->activePointY, ptrState->activePointZ);
        sprintf(msg2, "(np=%d, ne=%d, nf=%d)",
                ptrState->activePointI, ptrState->activePointJ, abs(ptrState->activePointK));
      }*/
      if (ptrState->mode == SCALARFIELD)
      {
        sprintf(msg, "%s=%g",
                zone->varnames[ptrState->scalarField], 
                ptrState->activePointF[ptrState->scalarField]);
      }
      else if (ptrState->mode == VECTORFIELD)
      {
        sprintf(msg, "%s=%g, %s=%g, %s=%g",
                zone->varnames[ptrState->vectorField1], 
                zone->f[ptrState->vectorField1][ptrState->activePointI],
                zone->varnames[ptrState->vectorField2],
                zone->f[ptrState->vectorField2][ptrState->activePointI],
                zone->varnames[ptrState->vectorField3],
                zone->f[ptrState->vectorField3][ptrState->activePointI]);
      }
      else
      {
        switch (zu->eltType)
        {
          case 0:
            sprintf(msg,"%s (NODE)", zone->zoneName); break;
          case 1:
            sprintf(msg,"%s (BAR)", zone->zoneName); break;
          case 2:
            sprintf(msg,"%s (TRI)", zone->zoneName); break;
          case 3:
            sprintf(msg,"%s (QUAD)", zone->zoneName); break;
          case 4:
            sprintf(msg,"%s (TETRA)", zone->zoneName); break;
          case 5:
            sprintf(msg,"%s (PENTA)", zone->zoneName); break;
          case 6:
            sprintf(msg,"%s (PYRA)", zone->zoneName); break;
          case 7:
            sprintf(msg,"%s (HEXA)", zone->zoneName); break;
          case 10:
            sprintf(msg,"%s (NGON)", zone->zoneName); break;
          default:
            sprintf(msg,"%s (UNKNOWN)", zone->zoneName); break;
        }
      }
    }

    // Patch by the CB
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    double winx,winy,winz;
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluProject(x, y, z, modelview, projection,
               viewport, &winx, &winy, &winz);
    double posX,posY,posZ;
    double shtx, shty, shtz;
    double shnx, shny, shnz;
    gluUnProject(winx+1., winy, winz, modelview, projection, viewport, 
                 &posX, &posY, &posZ);
    shtx = x-posX; shty = y-posY; shtz = z-posZ;
    gluUnProject(winx, winy+1., winz, modelview, projection, viewport, 
                 &posX, &posY, &posZ);
    shnx = x-posX; shny = y-posY; shnz = z-posZ;
    renderStringWithShadow(x+5.*r*(2.2*vx+dirx), y+5.*r*(2.2*vy+diry), 
                           z+5.*r*(2.2*vz+dirz), FONT3, msg,
                           1., 0., 0., 1., 
                           1., 1., 1., 0.7,
                           shtx, shty, shtz,
                           shnx, shny, shnz,
                           1.);
    if (msg2[0] != '\0')
      renderStringWithShadow(x+5*r*(2.2*vx+0.2*dirx), y+5*r*(2.2*vy+0.2*diry), 
                           z+5*r*(2.2*vz+0.2*dirz), FONT3, msg2,
                           1., 0., 0., 1., 
                           1., 1., 1., 0.7,
                           shtx, shty, shtz,
                           shnx, shny, shnz,
                           1.);
  }
 
  /* End of HUD */

  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_POLYGON_OFFSET_FILL);

  glLineWidth(1.);
  glColor4f(1., 1., 1., 1.);
}
