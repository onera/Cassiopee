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
#include "CompGeom/compGeom.h"

//==============================================================================
/*
  Retourne le max, le min, le ratio max sur min pour chaque edge
  d'un maillage
  IN: xt, yt, zt: coordonnees du maillage
  IN: im,jm,km: si maillage structure
  IN: cn: connectivite si maillage non structure
  IN: eltType: type d'elements
  IN: dim: 1/2/3
  IN: type=0 (max edge of cells), =1 (min edge of cells), =2 (ratio)
  OUT: out 
  Retourne 1 (OK) 0 (FAILED)
*/
//==============================================================================
E_Int K_COMPGEOM::getEdgeLength(E_Float* xt, E_Float* yt, E_Float* zt,
                                E_Int im, E_Int jm, E_Int km, 
                                K_FLD::FldArrayI* cn, char* eltType,
                                E_Int dim, E_Int type,
                                E_Float* out)
{
  E_Int indA, indB, indC, indD, indE, indF, indG, indH, indcell;
  E_Int ind1, ind2;
  E_Float dx, dy, dz, l;
  E_Float lout, lout2;

  if (im != -1) // structure
  {
    E_Int im1 = im-1;
    E_Int jm1 = jm-1;
    E_Int km1 = km-1;
    if (im == 1) im1 = 1;
    if (jm == 1) jm1 = 1;
    if (km == 1) km1 = 1;

    E_Int imjm = im*jm;
    E_Int imn, jmn, imc, imcjmc;
    if (dim == 2)
    {
      if (km1 == 1 && im1 > 1 && jm1 > 1) {imn = im; jmn = jm; imc = im1;}
      else if (im1 == 1 && jm1 > 1 && km1 > 1) {imn = jm; jmn = km; imc = jm1; }
      else if (jm1 == 1 && im1 > 1 && km1 > 1) {imn = km; jmn = im; imc = km1;}
      else 
      {
        printf("Error: getEdgeLength: constant direction for 2D pb cannot be determined.");
        return 0; // FAILED
      }

      if (type == 0) // max
      {
        for (E_Int j = 0; j < jmn-1; j++)
          for (E_Int i = 0; i < imn-1; i++)
          {
            indA = i+j*imn; indB = indA+1; indC = indB+imn; indD = indA+imn;
            indcell = i + j*imc;
            dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB 
            dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
            lout = l;
            
            dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD]; // CD
            dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_max(l,lout);

            dx = xt[indA]-xt[indD]; dy = yt[indA]-yt[indD]; // AD
            dz = zt[indA]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_max(l,lout);

            dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; // BC
            dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_max(l,lout);

            out[indcell] = sqrt(lout);
          }
      }
      else if (type == 1) // min
      {
        for (E_Int j = 0; j < jmn-1; j++)
          for (E_Int i = 0; i < imn-1; i++)
          {
            indA = i+j*imn; indB = indA+1; indC = indB+imn; indD = indA+imn;
            indcell = i + j*imc;
            dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB 
            dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
            lout = l;
            
            dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD]; // CD
            dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_min(l,lout);

            dx = xt[indA]-xt[indD]; dy = yt[indA]-yt[indD]; // AD
            dz = zt[indA]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_min(l,lout);

            dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; // BC
            dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_min(l,lout);

            out[indcell] = sqrt(lout);
          }
      }
      else if (type == 2) // ratio
      {
        for (E_Int j = 0; j < jmn-1; j++)
          for (E_Int i = 0; i < imn-1; i++)
          {
            indA = i+j*imn; indB = indA+1; indC = indB+imn; indD = indA+imn;
            indcell = i + j*imc;
            dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB 
            dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
            lout = l; lout2 = l;
            
            dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD]; // CD
            dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);

            dx = xt[indA]-xt[indD]; dy = yt[indA]-yt[indD]; // AD
            dz = zt[indA]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);

            dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; // BC
            dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
            lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);

            out[indcell] = sqrt(lout/lout2);
          }
      }
    }      
    else // dim3
    {
      E_Int edges1[12]; E_Int edges2[12];
      imcjmc = im1*jm1; imc = im1;
      if (type == 0)
      {
        for (E_Int k = 0; k < km-1; k++)
          for (E_Int j = 0; j < jm-1; j++)
            for (E_Int i = 0; i < im-1; i++)
            {
              indcell = i + j*imc + k*imcjmc;
              indA = i+j*im+k*imjm; indB = indA+1;
              indC = indB+im; indD = indA+im;
              indE = indA + imjm; indF = indB + imjm;
              indG = indC + imjm; indH = indD + imjm;        
              edges1[0] = indA; edges2[0] = indB; //AB
              edges1[1] = indD; edges2[1] = indC; //DC
              edges1[2] = indA; edges2[2] = indD; //AD
              edges1[3] = indB; edges2[3] = indC; //BC              
              edges1[4] = indE; edges2[4] = indF; //EF
              edges1[5] = indG; edges2[5] = indH; //GH
              edges1[6] = indE; edges2[6] = indH; //EH
              edges1[7] = indF; edges2[7] = indG; //FG
              edges1[8] = indA; edges2[8] = indE; //AE
              edges1[9] = indD; edges2[9] = indH; //DH
              edges1[10] = indB; edges2[10] = indF; //BF
              edges1[11] = indC; edges2[11] = indG; //CG
              
              lout =-K_CONST::E_MAX_FLOAT;
              for (E_Int ii = 0; ii < 12; ii++)
              {
                ind1 = edges1[ii]; ind2 = edges2[ii];
                dx = xt[ind1]-xt[ind2];
                dy = yt[ind1]-yt[ind2];
                dz = zt[ind1]-zt[ind2];
                l = dx*dx+dy*dy+dz*dz;
                lout = K_FUNC::E_max(l,lout);
              }
              out[indcell] = sqrt(lout);        
            }
      }
      else if (type == 1)
      {
        for (E_Int k = 0; k < km-1; k++)
          for (E_Int j = 0; j < jm-1; j++)
            for (E_Int i = 0; i < im-1; i++)
            {
              indcell = i + j*imc + k*imcjmc;
              indA = i+j*im+k*imjm; indB = indA+1;
              indC = indB+im; indD = indA+im;
              indE = indA + imjm; indF = indB + imjm;
              indG = indC + imjm; indH = indD + imjm;        
              edges1[0] = indA; edges2[0] = indB; //AB
              edges1[1] = indD; edges2[1] = indC; //DC
              edges1[2] = indA; edges2[2] = indD; //AD
              edges1[3] = indB; edges2[3] = indC; //BC              
              edges1[4] = indE; edges2[4] = indF; //EF
              edges1[5] = indG; edges2[5] = indH; //GH
              edges1[6] = indE; edges2[6] = indH; //EH
              edges1[7] = indF; edges2[7] = indG; //FG
              edges1[8] = indA; edges2[8] = indE; //AE
              edges1[9] = indD; edges2[9] = indH; //DH
              edges1[10] = indB; edges2[10] = indF; //BF
              edges1[11] = indC; edges2[11] = indG; //CG
              
              lout =+K_CONST::E_MAX_FLOAT;
              for (E_Int ii = 0; ii < 12; ii++)
              {
                ind1 = edges1[ii]; ind2 = edges2[ii];
                dx = xt[ind1]-xt[ind2];
                dy = yt[ind1]-yt[ind2];
                dz = zt[ind1]-zt[ind2];
                l = dx*dx+dy*dy+dz*dz;
                lout = K_FUNC::E_min(l,lout);
              }
              out[indcell] = sqrt(lout);        
            }
      }
      else if (type == 2)
      {
        for (E_Int k = 0; k < km-1; k++)
          for (E_Int j = 0; j < jm-1; j++)
            for (E_Int i = 0; i < im-1; i++)
            {
              indcell = i + j*imc + k*imcjmc;
              indA = i+j*im+k*imjm; indB = indA+1;
              indC = indB+im; indD = indA+im;
              indE = indA + imjm; indF = indB + imjm;
              indG = indC + imjm; indH = indD + imjm;        
              edges1[0] = indA; edges2[0] = indB; //AB
              edges1[1] = indD; edges2[1] = indC; //DC
              edges1[2] = indA; edges2[2] = indD; //AD
              edges1[3] = indB; edges2[3] = indC; //BC              
              edges1[4] = indE; edges2[4] = indF; //EF
              edges1[5] = indG; edges2[5] = indH; //GH
              edges1[6] = indE; edges2[6] = indH; //EH
              edges1[7] = indF; edges2[7] = indG; //FG
              edges1[8] = indA; edges2[8] = indE; //AE
              edges1[9] = indD; edges2[9] = indH; //DH
              edges1[10] = indB; edges2[10] = indF; //BF
              edges1[11] = indC; edges2[11] = indG; //CG
              
              lout =-K_CONST::E_MAX_FLOAT;
              lout2 =+K_CONST::E_MAX_FLOAT;
              for (E_Int ii = 0; ii < 12; ii++)
              {
                ind1 = edges1[ii]; ind2 = edges2[ii];
                dx = xt[ind1]-xt[ind2];
                dy = yt[ind1]-yt[ind2];
                dz = zt[ind1]-zt[ind2];
                l = dx*dx+dy*dy+dz*dz;
                lout = K_FUNC::E_max(l,lout);
                lout2 = K_FUNC::E_min(l,lout2);
              }
              out[indcell] = sqrt(lout/lout2);        
            }
      }
    }
    return 1; // OK
  }

  E_Int nedges = 0;
  if (strcmp(eltType, "BAR") == 0 || strcmp(eltType, "NODE") == 0) 
  {
    printf("Error: getEdgeLength: not valid for NODE or BAR elements.");
    return 0;    
  }
  else if (strcmp(eltType, "TRI") == 0) nedges = 3; 
  else if (strcmp(eltType, "QUAD") == 0 ) nedges = 4;
  else if (strcmp(eltType, "TETRA") == 0) nedges = 6;
  else if (strcmp(eltType, "HEXA") == 0) nedges = 12;
  else if (strcmp(eltType, "PENTA") == 0) nedges = 9;
  else if (strcmp(eltType, "PYRA") == 0) nedges = 8;
  else if (strcmp(eltType, "NGON") == 0) nedges = -1;
  else
  {
    printf("Error: getEdgeLength: unknown type of element.");
    return 0;
  }

  E_Int nelts = cn->getSize();

  // CAS TRI
  if (nedges == 3) 
  {
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3);
    if (type == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1;
        
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l;
        
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; //BC
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
      
        dx = xt[indC]-xt[indA]; dy = yt[indC]-yt[indA]; // CA
        dz = zt[indC]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        out[et] = sqrt(lout);
      }
    }
    else if (type == 1)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1;
        
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l;
        
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; //BC
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
      
        dx = xt[indC]-xt[indA]; dy = yt[indC]-yt[indA]; // CA
        dz = zt[indC]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        out[et] = sqrt(lout);
      }
    }
    else if (type == 2)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1;
        
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l; lout2 = l;
        
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; //BC
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);
      
        dx = xt[indC]-xt[indA]; dy = yt[indC]-yt[indA]; // CA
        dz = zt[indC]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);
        out[et] = sqrt(lout/lout2);
      }
    }
  }
  else if (nedges == 4) // QUAD
  {
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3);
    E_Int* cn4 = cn->begin(4);
    if (type == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l;
        
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; //BC
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        
        dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD]; // CD
        dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        
        dx = xt[indD]-xt[indA]; dy = yt[indD]-yt[indA]; //AD
        dz = zt[indD]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        out[et] = sqrt(lout);
      }
    }
    else if (type == 1)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l;
        
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; //BC
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        
        dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD]; // CD
        dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        
        dx = xt[indD]-xt[indA]; dy = yt[indD]-yt[indA]; //AD
        dz = zt[indD]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        out[et] = sqrt(lout);
      }
    }
    else if (type == 2)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA]; // AB
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l; lout2 = l;
        
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB]; //BC
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);
        
        dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD]; // CD
        dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);
        
        dx = xt[indD]-xt[indA]; dy = yt[indD]-yt[indA]; //AD
        dz = zt[indD]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout); lout2 = K_FUNC::E_min(l,lout2);
        out[et] = sqrt(lout/lout2);
      }
    }
  }
  else if (nedges == 6) //TETRA
  {
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3);
    E_Int* cn4 = cn->begin(4);
    if (type == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        // AB
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA];
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l;
        //BC
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB];
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        // AC
        dx = xt[indC]-xt[indA]; dy = yt[indC]-yt[indA];
        dz = zt[indC]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        //AD
        dx = xt[indD]-xt[indA]; dy = yt[indD]-yt[indA];
        dz = zt[indD]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        // CD
        dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD];
        dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        // BD
        dx = xt[indB]-xt[indD]; dy = yt[indB]-yt[indD];
        dz = zt[indB]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        out[et] = sqrt(lout);
      }
    }
    else if (type == 1)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        // AB
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA];
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l;
        //BC
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB];
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        // AC
        dx = xt[indC]-xt[indA]; dy = yt[indC]-yt[indA];
        dz = zt[indC]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        //AD
        dx = xt[indD]-xt[indA]; dy = yt[indD]-yt[indA];
        dz = zt[indD]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        // CD
        dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD];
        dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        // BD
        dx = xt[indB]-xt[indD]; dy = yt[indB]-yt[indD];
        dz = zt[indB]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        out[et] = sqrt(lout);
      }
    }
    else if (type == 2)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        // AB
        dx = xt[indB]-xt[indA]; dy = yt[indB]-yt[indA];
        dz = zt[indB]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = l; lout2 = l;
        //BC
        dx = xt[indC]-xt[indB]; dy = yt[indC]-yt[indB];
        dz = zt[indC]-zt[indB]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        lout2 = K_FUNC::E_min(l,lout2);
        // AC
        dx = xt[indC]-xt[indA]; dy = yt[indC]-yt[indA];
        dz = zt[indC]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        lout2 = K_FUNC::E_min(l,lout2);
        //AD
        dx = xt[indD]-xt[indA]; dy = yt[indD]-yt[indA];
        dz = zt[indD]-zt[indA]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        lout2 = K_FUNC::E_min(l,lout2);
        // CD
        dx = xt[indC]-xt[indD]; dy = yt[indC]-yt[indD];
        dz = zt[indC]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        lout2 = K_FUNC::E_min(l,lout2);
        // BD
        dx = xt[indB]-xt[indD]; dy = yt[indB]-yt[indD];
        dz = zt[indB]-zt[indD]; l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        lout2 = K_FUNC::E_min(l,lout2);
        out[et] = sqrt(lout/lout2);
      }
    }
  }

  else if (nedges == 12) // HEXA        
  {
    E_Int edges1[12];
    E_Int edges2[12];
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3);
    E_Int* cn4 = cn->begin(4);
    E_Int* cn5 = cn->begin(5);
    E_Int* cn6 = cn->begin(6);
    E_Int* cn7 = cn->begin(7);
    E_Int* cn8 = cn->begin(8);
    if (type == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        indE = cn5[et]-1; indF = cn6[et]-1; indG = cn7[et]-1; indH = cn8[et]-1;
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indD; edges2[1] = indC; //DC
        edges1[2] = indA; edges2[2] = indD; //AD
        edges1[3] = indB; edges2[3] = indC; //BC              
        edges1[4] = indE; edges2[4] = indF; //EF
        edges1[5] = indG; edges2[5] = indH; //GH
        edges1[6] = indE; edges2[6] = indH; //EH
        edges1[7] = indF; edges2[7] = indG; //FG
        edges1[8] = indA; edges2[8] = indE; //AE
        edges1[9] = indD; edges2[9] = indH; //DH
        edges1[10] = indB; edges2[10] = indF; //BF
        edges1[11] = indC; edges2[11] = indG; //CG
        
        lout =-K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2]; dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2]; l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
        }
        // Edge Ratio 
        out[et] = sqrt(lout);    
      }
    }
    else if (type == 1)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        indE = cn5[et]-1; indF = cn6[et]-1; indG = cn7[et]-1; indH = cn8[et]-1;
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indD; edges2[1] = indC; //DC
        edges1[2] = indA; edges2[2] = indD; //AD
        edges1[3] = indB; edges2[3] = indC; //BC              
        edges1[4] = indE; edges2[4] = indF; //EF
        edges1[5] = indG; edges2[5] = indH; //GH
        edges1[6] = indE; edges2[6] = indH; //EH
        edges1[7] = indF; edges2[7] = indG; //FG
        edges1[8] = indA; edges2[8] = indE; //AE
        edges1[9] = indD; edges2[9] = indH; //DH
        edges1[10] = indB; edges2[10] = indF; //BF
        edges1[11] = indC; edges2[11] = indG; //CG
        
        lout = K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2]; dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2]; l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_min(l,lout);
        }
        // Edge Ratio 
        out[et] = sqrt(lout);    
      }
    }
    else if (type == 2)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
        indE = cn5[et]-1; indF = cn6[et]-1; indG = cn7[et]-1; indH = cn8[et]-1;
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indD; edges2[1] = indC; //DC
        edges1[2] = indA; edges2[2] = indD; //AD
        edges1[3] = indB; edges2[3] = indC; //BC              
        edges1[4] = indE; edges2[4] = indF; //EF
        edges1[5] = indG; edges2[5] = indH; //GH
        edges1[6] = indE; edges2[6] = indH; //EH
        edges1[7] = indF; edges2[7] = indG; //FG
        edges1[8] = indA; edges2[8] = indE; //AE
        edges1[9] = indD; edges2[9] = indH; //DH
        edges1[10] = indB; edges2[10] = indF; //BF
        edges1[11] = indC; edges2[11] = indG; //CG
        
        lout =-K_CONST::E_MAX_FLOAT;
        lout2 = K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2]; dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2]; l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
          lout2 = K_FUNC::E_min(l,lout2);
        }
        // Edge Ratio 
        out[et] = sqrt(lout/lout2);    
      }
    }
  }
  else if (nedges == 9) //PENTA
  {
    E_Int edges1[9];
    E_Int edges2[9];
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3);
    E_Int* cn4 = cn->begin(4);
    E_Int* cn5 = cn->begin(5);
    E_Int* cn6 = cn->begin(6);
    
    if (type == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; 
        indD = cn4[et]-1; indE = cn5[et]-1; indF = cn6[et]-1; 
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indB; edges2[1] = indC; //BC
        edges1[2] = indC; edges2[2] = indA; //CA
        edges1[3] = indD; edges2[3] = indE; //DE              
        edges1[4] = indE; edges2[4] = indF; //EF
        edges1[5] = indF; edges2[5] = indD; //FD
        edges1[6] = indA; edges2[6] = indD; //AD
        edges1[7] = indC; edges2[7] = indF; //CF
        edges1[8] = indB; edges2[8] = indE; //BE
        
        lout =-K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2];
          dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
        } 
        out[et] = sqrt(lout);
      }
    }
    else if (type == 1)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; 
        indD = cn4[et]-1; indE = cn5[et]-1; indF = cn6[et]-1; 
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indB; edges2[1] = indC; //BC
        edges1[2] = indC; edges2[2] = indA; //CA
        edges1[3] = indD; edges2[3] = indE; //DE              
        edges1[4] = indE; edges2[4] = indF; //EF
        edges1[5] = indF; edges2[5] = indD; //FD
        edges1[6] = indA; edges2[6] = indD; //AD
        edges1[7] = indC; edges2[7] = indF; //CF
        edges1[8] = indB; edges2[8] = indE; //BE
        
        lout = K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2];
          dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_min(l,lout);
        } 
        out[et] = sqrt(lout);
      }
    }
    else if (type == 2)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; 
        indD = cn4[et]-1; indE = cn5[et]-1; indF = cn6[et]-1; 
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indB; edges2[1] = indC; //BC
        edges1[2] = indC; edges2[2] = indA; //CA
        edges1[3] = indD; edges2[3] = indE; //DE              
        edges1[4] = indE; edges2[4] = indF; //EF
        edges1[5] = indF; edges2[5] = indD; //FD
        edges1[6] = indA; edges2[6] = indD; //AD
        edges1[7] = indC; edges2[7] = indF; //CF
        edges1[8] = indB; edges2[8] = indE; //BE
        
        lout =-K_CONST::E_MAX_FLOAT;
        lout2 = K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2];
          dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
          lout2 = K_FUNC::E_min(l,lout2);
        } 
        out[et] = sqrt(lout/lout2);
      }
    }
  }
  
  else if (nedges == 8) //PYRA
  {
    E_Int edges1[8];
    E_Int edges2[8];
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3);
    E_Int* cn4 = cn->begin(4);
    E_Int* cn5 = cn->begin(5);
    
    if (type == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; 
        indD = cn4[et]-1; indE = cn5[et]-1;
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indB; edges2[1] = indC; //BC
        edges1[2] = indC; edges2[2] = indD; //CD
        edges1[3] = indD; edges2[3] = indA; //DA              
        edges1[4] = indA; edges2[4] = indE; //AE
        edges1[5] = indB; edges2[5] = indE; //BE
        edges1[6] = indC; edges2[6] = indE; //CE
        edges1[7] = indD; edges2[7] = indE; //DE
        
        lout =-K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2];
          dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
        }
        out[et] = sqrt(lout);    
      }
    }
    else if (type == 1)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; 
        indD = cn4[et]-1; indE = cn5[et]-1;
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indB; edges2[1] = indC; //BC
        edges1[2] = indC; edges2[2] = indD; //CD
        edges1[3] = indD; edges2[3] = indA; //DA              
        edges1[4] = indA; edges2[4] = indE; //AE
        edges1[5] = indB; edges2[5] = indE; //BE
        edges1[6] = indC; edges2[6] = indE; //CE
        edges1[7] = indD; edges2[7] = indE; //DE
        
        lout = K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2];
          dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_min(l,lout);
        } 
        out[et] = sqrt(lout);    
      }
    }
    else if (type == 2)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; 
        indD = cn4[et]-1; indE = cn5[et]-1;
        edges1[0] = indA; edges2[0] = indB; //AB
        edges1[1] = indB; edges2[1] = indC; //BC
        edges1[2] = indC; edges2[2] = indD; //CD
        edges1[3] = indD; edges2[3] = indA; //DA              
        edges1[4] = indA; edges2[4] = indE; //AE
        edges1[5] = indB; edges2[5] = indE; //BE
        edges1[6] = indC; edges2[6] = indE; //CE
        edges1[7] = indD; edges2[7] = indE; //DE
        
        lout =-K_CONST::E_MAX_FLOAT;
        lout2 = K_CONST::E_MAX_FLOAT;
        for (E_Int ii = 0; ii < nedges; ii++)
        {
          ind1 = edges1[ii]; ind2 = edges2[ii];
          dx = xt[ind1]-xt[ind2];
          dy = yt[ind1]-yt[ind2];
          dz = zt[ind1]-zt[ind2];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
          lout2 = K_FUNC::E_min(l,lout2);
        } 
        out[et] = sqrt(lout/lout2);
      }
    }
  }

  else if (nedges == -1) // NGON
  {
    E_Int* cnp = cn->begin(); // pointeur sur la connectivite NGon
    E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
    E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
    E_Int nfaces = cnp[0];
    E_Int nvert, noface;

    // for all the faces computes the min and max length of edges
    E_Int* cFV = cnp+2;

    if (type == 0)
    {
      E_Float* loutp = new E_Float [nfaces];
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        nvert = cFV[0];
        lout =-K_CONST::E_MAX_FLOAT;
        for (E_Int i = 1; i < nvert; i++)
        {
          indA = cFV[i]-1; indB = cFV[i+1]-1;
          dx = xt[indA]-xt[indB];
          dy = yt[indA]-yt[indB];
          dz = zt[indA]-zt[indB];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
        }
        // on boucle
        indA = cFV[nvert]-1; indB = cFV[1]-1;
        dx = xt[indA]-xt[indB];
        dy = yt[indA]-yt[indB];
        dz = zt[indA]-zt[indB];
        l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        
        loutp[nof] = lout;
        cFV += nvert+1;
      }
      E_Int* ptr = cnp+sizeFN+4;
      
      // for any element, ratio of max/min edges for all its faces
      for (E_Int noe = 0; noe < nelts; noe++)
      {
        nfaces = ptr[0];
        lout = K_CONST::E_MAX_FLOAT;
        for (E_Int j = 1; j <= nfaces; j++)
        {
          noface = ptr[j]-1;
          lout = K_FUNC::E_min(lout,loutp[noface]);
        }
        ptr += nfaces+1;
        out[noe] = sqrt(lout);
      }
      delete [] loutp;
    }
    else if (type == 1)
    {
      E_Float* loutp = new E_Float [nfaces];
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        E_Int nvert = cFV[0];
        lout = K_CONST::E_MAX_FLOAT;
        for (E_Int i = 1; i < nvert; i++)
        {
          indA = cFV[i]-1; indB = cFV[i+1]-1;
          dx = xt[indA]-xt[indB];
          dy = yt[indA]-yt[indB];
          dz = zt[indA]-zt[indB];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_min(l,lout);
        }
        // on boucle
        indA = cFV[nvert]-1; indB = cFV[1]-1;
        dx = xt[indA]-xt[indB];
        dy = yt[indA]-yt[indB];
        dz = zt[indA]-zt[indB];
        l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_min(l,lout);
        
        loutp[nof] = lout;
        cFV += nvert+1;
      }
      E_Int* ptr = cnp+sizeFN+4;
      
      // for any element, ratio of max/min edges for all its faces
      for (E_Int noe = 0; noe < nelts; noe++)
      {
        nfaces = ptr[0];
        lout = K_CONST::E_MAX_FLOAT;
        for (E_Int j = 1; j <= nfaces; j++)
        {
          noface = ptr[j]-1;
          lout = K_FUNC::E_min(lout,loutp[noface]);
        }
        ptr += nfaces+1;
        out[noe] = sqrt(lout);
      }
      delete [] loutp;
    }
    else if (type == 2)
    {
      E_Float* loutp = new E_Float [nfaces];
      E_Float* lout2p = new E_Float [nfaces];
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        E_Int nvert = cFV[0];
        lout =-K_CONST::E_MAX_FLOAT;
        lout2 = K_CONST::E_MAX_FLOAT;
        for (E_Int i = 1; i < nvert; i++)
        {
          indA = cFV[i]-1; indB = cFV[i+1]-1;
          dx = xt[indA]-xt[indB];
          dy = yt[indA]-yt[indB];
          dz = zt[indA]-zt[indB];
          l = dx*dx+dy*dy+dz*dz;
          lout = K_FUNC::E_max(l,lout);
          lout2 = K_FUNC::E_min(l,lout2);
        }
        // on boucle
        indA = cFV[nvert]-1; indB = cFV[1]-1;
        dx = xt[indA]-xt[indB];
        dy = yt[indA]-yt[indB];
        dz = zt[indA]-zt[indB];
        l = dx*dx+dy*dy+dz*dz;
        lout = K_FUNC::E_max(l,lout);
        lout2 = K_FUNC::E_min(l,lout2);
        
        loutp[nof] = lout; lout2p[nof] = lout2;
        cFV += nvert+1;
      }
      E_Int* ptr = cnp+sizeFN+4;
      
      // for any element, ratio of max/min edges for all its faces
      for (E_Int noe = 0; noe < nelts; noe++)
      {
        nfaces = ptr[0];
        lout =-K_CONST::E_MAX_FLOAT;
        lout2 = K_CONST::E_MAX_FLOAT;
        for (E_Int j = 1; j <= nfaces; j++)
        {
          noface = ptr[j]-1;
          lout = K_FUNC::E_max(lout,loutp[noface]);
          lout2 = K_FUNC::E_min(lout2,lout2p[noface]);
        }
        ptr += nfaces+1;
        out[noe] = sqrt(lout/lout2);
      }
      delete [] loutp; delete [] lout2p;
    }

  }
    return 1;
}
