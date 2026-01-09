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

  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;

  for (plane = 0; plane < 6; plane++)
  {    
    nis = 0; nie = ni;
    njs = 0; nje = nj;
    nks = 0; nke = nk;
    
    switch (plane)
    {
      case 0:
        nis = zonep->iPlane;
        if (nis == -1)
        { nis = 0; nie = 1; }
        else if (nis == -2)
        { nis = 0; nie = 0; }
        else nie = nis+1;
        break;
            
      case 1:
        njs = zonep->jPlane;
        if (njs == -1)
        { njs = 0; nje = 1;}
        else if (njs == -2)
        { njs = 0; nje = 0; }
        else nje = njs+1;
        break;
          
      case 2:
        nks = zonep->kPlane;
        if (nks == -1)
        { nks = 0; nke = 1; }
        else if (nks == -2)
        { nks = 0; nke = 0; }
        else nke = nks+1;
        break;
        
      case 3:
        nis = zonep->iPlane;
        if (nis == -1)
        { nis = ni-1; nie = ni; }
        else
        { nie = 0; nje = 0; nke = 0;}
        break;
        
      case 4:
        njs = zonep->jPlane;
        if (njs == -1)
        { njs = nj-1; nje = nj; }
        else
        { nie = 0; nje = 0; nke = 0;}  
        break;
        
      case 5:
        nks = zonep->kPlane;
        if (nks == -1)
        { nks = nk-1; nke = nk; }
        else
        { nie = 0; nje = 0; nke = 0;}  
        break;
    }

    if (zonep->blank == 0) // No blanking
    {
      for (k = nks; k < nke; k += stepk)
        for (j = njs; j < nje; j += stepj)
        {
          glBegin(GL_LINE_STRIP);
          n1 = j*ni+k*nij;
          glVertex3d(x[n1+nis], y[n1+nis], z[n1+nis]);
          for (i = nis+stepi; i < nie; i = i+stepi)
          {
            n2 = n1+i;       
            glVertex3d(x[n2], y[n2], z[n2]);
          }
          glEnd();
        }

      for (k = nks; k < nke; k += stepk)
        for (i = nis; i < nie; i += stepi)
        {
          glBegin(GL_LINE_STRIP);
          n1 = i+k*nij;
          glVertex3d(x[n1+njs*ni], y[n1+njs*ni], z[n1+njs*ni]);
          for (j = njs+stepj; j < nje; j = j+stepj)
          {
            n2 = n1+j*ni;
            glVertex3d(x[n2], y[n2], z[n2]);
          }
          glEnd();
        }

      for (j = njs; j < nje; j += stepj)
        for (i = nis; i < nie; i += stepi)
        {
          glBegin(GL_LINE_STRIP);
          n1 = i+j*ni;
          glVertex3d(x[n1+nks*nij], y[n1+nks*nij], z[n1+nks*nij]);
          for (k = nks+stepk; k < nke; k = k+stepk)
          {
            n2 = n1+k*nij;
            glVertex3d(x[n2], y[n2], z[n2]);
          }
          glEnd();
        }
    }
    else // With blanking
    {
      glBegin(GL_LINES);

      for (k = nks; k < nke; k += stepk)
        for (j = njs; j < nje; j += stepj)
        {
          n1 = nis+j*ni+k*nij;
          ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
          for (i = nis; i < nie-stepi; i = i+stepi)
          {
            n2 = n1+stepi;
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              glVertex3d(x[n1], y[n1], z[n1]);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
            ret1 = ret2; n1 = n2;
          }
        }
      for (k = nks; k < nke; k += stepk)
        for (i = nis; i < nie; i += stepi)
        {
          n1 = i+njs*ni+k*nij;
          ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
          for (j = njs; j < nje-stepj; j = j+stepj)
          {
            n2 = n1+ni*stepj;
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              glVertex3d(x[n1], y[n1], z[n1]);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
            ret1 = ret2; n1 = n2;
          }
        }

      for (j = njs; j < nje; j += stepj)
        for (i = nis; i < nie; i += stepi)
        {
          n1 = i+j*ni+nks*nij;
          ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
          for (k = nks; k < nke-stepk; k = k+stepk)
          {
            n2 = n1+nij*stepk;
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              glVertex3d(x[n1], y[n1], z[n1]);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
            ret1 = ret2; n1 = n2;
          }
        }

      glEnd();
    }  
  } // planes
