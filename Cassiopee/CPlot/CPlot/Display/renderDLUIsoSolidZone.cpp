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
#include "../DataDL.h"

//=============================================================================
/*
  Display une zone en iso solides pour un champ scalaire.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield: le no du champ a afficher.
*/
//=============================================================================
void DataDL::renderUIsoSolidZone(UnstructZone* zonep, E_Int zonet, E_Int nofield)
{
  float offb;
  double blend;
  
  ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  
  // Blending
  blend = 1.;
#include "selection2.h"

  E_Int eltType0 = zonep->eltType[0];
  bool is1D = false;
  if ((eltType0 == 1) || (eltType0 == 10 && zonep->nelts1D > 0)) is1D = true;
  
#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend

  if (is1D)
  {
    if (curr != 0) _shaders[curr]->setUniform("lightOn", (int)0); // impose isoLight off on 1D meshes
  }
#endif
  
  glCallList(zoneImpl->_DLiso);

#ifdef __SHADERS__
  if (is1D)
  {
    if (ptrState->isoLight == 1 && ptrState->dim == 3)
    {
      if (curr != 0) _shaders[curr]->setUniform("lightOn", (int)1); // put back the isoLight value found in the CPlot state
    }
  }
#endif
}

//=============================================================================
/*
  Display une zone en iso solides pour un champ vectoriel.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield1, nofield2, nofield3: les no des champs a afficher.
*/
//=============================================================================
void DataDL::renderUIsoSolidZone(UnstructZone* zonep, E_Int zonet,
				 E_Int nofield1, E_Int nofield2, E_Int nofield3)
{
  E_Int i, n1, n2;
  float offb;
  double blend;
  E_Int ret1, ret2;

  ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  
  // Blending
  blend = 1.;
#include "selection2.h"

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

  double* x = zonep->x; double* y = zonep->y; double* z = zonep->z;
  E_Int ne = zonep->nec[0];
  E_Int eltType0 = zonep->eltType[0];
  E_Int* connect = zonep->connect[0];
  
  glCallList(zoneImpl->_DLiso);

  // Pour les BAR
  if (eltType0 == 1)
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
  if (eltType0 == 10 && zonep->nelts1D > 0)
  {
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        E_Int elt = zonep->posElts1D[i];
        E_Int* ptrelt = &connect[elt];
        E_Int face1 = ptrelt[1]-1;
        E_Int face2 = ptrelt[2]-1;
        E_Int posface1 = zonep->posFaces[face1];
        E_Int posface2 = zonep->posFaces[face2];
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
        E_Int elt = zonep->posElts1D[i];
        E_Int* ptrelt = &connect[elt];
        E_Int face1 = ptrelt[1]-1;
        E_Int face2 = ptrelt[2]-1;
        E_Int posface1 = zonep->posFaces[face1];
        E_Int posface2 = zonep->posFaces[face2];
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
