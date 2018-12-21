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
#define PLOT glVertex3d(x[n1], y[n1], z[n1]); \
  glVertex3d(x[n2], y[n2], z[n2]);
#define PLOTB ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);             \
  if (ret1*ret2 != 0) { PLOT; }

  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;
  int* connect = zonep->connect;
  
  // Grid dimensions
  int ne = zonep->ne;
  int ne2 = 2*ne;
  int ne3 = 3*ne;
  int eltType = zonep->eltType;
  int nd, l;

  glBegin(GL_LINES);
  
  if (zonep->blank == 0)
  {
    switch (eltType)
    {
      case 1: // BAR
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOT;
        }
        break;
        
      case 2: // TRI
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;       
          PLOT;
          
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOT;
          
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          PLOT;
        }
        break;
        
      case 3: // QUAD
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOT;
          
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOT;
          
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne3]-1;
          PLOT;
          
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOT;
        }
        break;
        
      case 4: // TETRA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOT;
          
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOT;
          
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          PLOT;
          
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOT;
          
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne]-1;
          PLOT;
          
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne2]-1;
          PLOT;
        }
        break;

      case 5: // PENTA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOT;
          
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOT;
          
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          PLOT;
          
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
          
          n1 = connect[i+ne*4]-1;
          n2 = connect[i+ne*5]-1;
          PLOT;
          
          n1 = connect[i+ne*5]-1;
          n2 = connect[i+ne3]-1;
          PLOT;
          
          n1 = connect[i]-1;
          n2 = connect[i+ne3]-1;
          PLOT;
          
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
          
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne*5]-1;
          PLOT;
        }
        break;

      case 6: // PYRA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOT;              
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOT;              
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne3]-1;
          PLOT;
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOT;              
          n1 = connect[i]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;           
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
        }
        break;

      case 7: // HEXA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOT;              
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOT;              
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne3]-1;
          PLOT;
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOT;              
          n1 = connect[i+ne*4]-1;
          n2 = connect[i+ne*5]-1;
          PLOT;           
          n1 = connect[i+ne*5]-1;
          n2 = connect[i+ne*6]-1;
          PLOT;
          n1 = connect[i+ne*6]-1;
          n2 = connect[i+ne*7]-1;
          PLOT;
          n1 = connect[i+ne*7]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
          n1 = connect[i]-1;
          n2 = connect[i+ne*4]-1;
          PLOT;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne*5]-1;
          PLOT;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne*6]-1;
          PLOT;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne*7]-1;
          PLOT;
        }
        break;

      case 10: // NGON
      {
        int nf = NFACES(connect); // nbre de faces
        int *ptrface = PTRFACES(connect);
        for (i = 0; i < nf; i++)
        {
          nd = ptrface[0]; // nbre de noeuds de la face
          for (l = 0; l < nd-1; l++)
          {
            n1 = ptrface[l+1]-1; n2 = ptrface[l+2]-1;
            PLOT;
          }
          n1 = ptrface[nd]-1; n2 = ptrface[1]-1;
          PLOT;
          ptrface += nd+1;
        }

        // Elements 1D
        int elt;
        for (i = 0; i < zonep->nelts1D; i++)
        {
          elt = zonep->posElts1D[i];
          int* ptrelt = &connect[elt];
          int face = ptrelt[1]-1; // indice de la face
          int* ptrface = &connect[zonep->posFaces[face]];
          n1 = ptrface[1]-1;
          face = ptrelt[2]-1; // indice de la face
          ptrface = &connect[zonep->posFaces[face]];
          n2 = ptrface[1]-1;
          PLOT;
        }
      }
      break;
    }
  }
  else // With blanking
  {
    switch (eltType)
    {
      case 0: // NODE
        // done later
        break;
      case 1: // BAR
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
        }
        break;
        
      case 2: // TRI
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          PLOTB;
        }
        break;
        
      case 3: // QUAD
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne3]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOTB;
        }
        break;
        
      case 4: // TETRA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;
        }
        break;
        
      case 5: // PENTA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
          n1 = connect[i+ne*4]-1;
          n2 = connect[i+ne*5]-1;
          PLOTB;
          n1 = connect[i+ne*5]-1;
          n2 = connect[i+ne3]-1;
          PLOTB;
          n1 = connect[i]-1;
          n2 = connect[i+ne3]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne*5]-1;
          PLOTB;
        }
        break;

      case 6: // PYRA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;              
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;              
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne3]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOTB;              
          n1 = connect[i]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;           
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
        }
        break;

      case 7: // HEXA
        for (i = 0; i < ne; i++)
        {
          n1 = connect[i]-1;
          n2 = connect[i+ne]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne2]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne3]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i]-1;
          PLOTB;
          n1 = connect[i+ne*4]-1;
          n2 = connect[i+ne*5]-1;
          PLOTB;
          n1 = connect[i+ne*5]-1;
          n2 = connect[i+ne*6]-1;
          PLOTB;
          n1 = connect[i+ne*6]-1;
          n2 = connect[i+ne*7]-1;
          PLOTB;
          n1 = connect[i+ne*7]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
          n1 = connect[i]-1;
          n2 = connect[i+ne*4]-1;
          PLOTB;
          n1 = connect[i+ne]-1;
          n2 = connect[i+ne*5]-1;
          PLOTB;
          n1 = connect[i+ne2]-1;
          n2 = connect[i+ne*6]-1;
          PLOTB;
          n1 = connect[i+ne3]-1;
          n2 = connect[i+ne*7]-1;
          PLOTB;
        }
        break;

      case 10: // NGON
      {
        int nf = NFACES(connect); // nbre de faces
        int *ptrface = PTRFACES(connect);
        for (i = 0; i < nf; i++)
        {
          nd = ptrface[0]; // nbre de noeuds de la face
          for (l = 0; l < nd-1; l++)
          {
            n1 = ptrface[l+1]-1; n2 = ptrface[l+2]-1;
            PLOTB;
          }
          n1 = ptrface[nd]-1; n2 = ptrface[1]-1;
          PLOTB;
          ptrface += nd+1;
        }

        // Elements 1D
        int elt;
        for (i = 0; i < zonep->nelts1D; i++)
        {
          elt = zonep->posElts1D[i];
          int* ptrelt = &connect[elt];
          int face = ptrelt[1]-1; // indice de la face
          int* ptrface = &connect[zonep->posFaces[face]];
          n1 = ptrface[1]-1;
          face = ptrelt[2]-1; // indice de la face
          ptrface = &connect[zonep->posFaces[face]];
          n2 = ptrface[1]-1;
          PLOTB;
        }
      }
      break;
    }
  }
  glEnd();
  //glLineWidth(1.);

