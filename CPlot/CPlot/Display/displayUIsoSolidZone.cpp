/*    
    Copyright 2013-2018 Onera.

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
  Display une zone en iso solides pour les champs scalaires.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield: le no du champ a afficher.
*/
//=============================================================================
void Data::displayUIsoSolidZone(UnstructZone* zonep, int zonet,
                                int nofield)
{
  int i, n1, n2, n3, n4;
  float r, g, b, offb;
  double blend;
  int ret1, ret2, ret3, ret4, ff;

  // Blending
  blend = 1.;
#include "selection2.h"
 
#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

#include "displayUIsoSolidZone.h"

  // Pour les BARS
  if (zonep->eltType == 1)
  {
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
        glColor4f(0., 0., 0.+offb, blend);   
        glVertex3d(x[n1], y[n1], z[n1]);
        getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);
        glColor4f(0., 0., 0.+offb, blend);   
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        
        if (ret1*ret2 != 0)
        {
          getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n1], y[n1], z[n1]);
          getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
  }

  // Pour les NGONS 1D
  if (zonep->eltType == 10 && zonep->nelts1D > 0)
  {
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        int elt = zonep->posElts1D[i];
        int* ptrelt = &connect[elt];
        int face1 = ptrelt[1]-1;
        int face2 = ptrelt[2]-1;
        int posface1 = zonep->posFaces[face1];
        int posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        glColor4f(0., 0., 0.+offb, blend); 
        glVertex3d(x[n1], y[n1], z[n1]);
        glColor4f(0., 0., 0.+offb, blend); 
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        int elt = zonep->posElts1D[i];
        int* ptrelt = &connect[elt];
        int face1 = ptrelt[1]-1;
        int face2 = ptrelt[2]-1;
        int posface1 = zonep->posFaces[face1];
        int posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        if (ret1*ret2 != 0)
        {
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n1], y[n1], z[n1]);
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
  }
}

//=============================================================================
/*
  Display une zone en iso solides pour les champs vectoriels.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield1, nofield2, nofield3: les no des champs a afficher.
*/
//=============================================================================
void Data::displayUIsoSolidZone(UnstructZone* zonep, int zonet,
                                int nofield1, int nofield2, int nofield3)
{
  int i, n1, n2, n3, n4;
  float r, g, b, offb;
  double blend;
  int ret1, ret2, ret3, ret4, ff;

  // Blending
  blend = 1.;
#include "selection2.h"

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

#undef PLOTTRI
#undef PLOTTRI2
#undef PLOTQUAD
#undef PLOTQUAD2
#undef PLOTNGON
#include "displayUVectSolidZone.h"

  // Pour les BARS
  if (zonep->eltType == 1)
  {
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        glColor4f(0., 0., 0.+offb, blend);   
        glVertex3d(x[n1], y[n1], z[n1]);
        glColor4f(0., 0., 0.+offb, blend);   
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        
        if (ret1*ret2 != 0)
        {
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n1], y[n1], z[n1]);
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
  }

  // Pour les NGONS 1D
  if (zonep->eltType == 10 && zonep->nelts1D > 0)
  {
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        int elt = zonep->posElts1D[i];
        int* ptrelt = &connect[elt];
        int face1 = ptrelt[1]-1;
        int face2 = ptrelt[2]-1;
        int posface1 = zonep->posFaces[face1];
        int posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        glColor4f(0., 0., 0.+offb, blend); 
        glVertex3d(x[n1], y[n1], z[n1]);
        glColor4f(0., 0., 0.+offb, blend); 
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        int elt = zonep->posElts1D[i];
        int* ptrelt = &connect[elt];
        int face1 = ptrelt[1]-1;
        int face2 = ptrelt[2]-1;
        int posface1 = zonep->posFaces[face1];
        int posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        if (ret1*ret2 != 0)
        {
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n1], y[n1], z[n1]);
          glColor4f(0., 0., 0.+offb, blend); 
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
  }
}
