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

// Formated xfig file support

# include <stdio.h>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include <vector>
# include "Def/DefFunction.h"
# include "Connect/connect.h"
# include "CompGeom/compGeom.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* xfigread 
   Read xfig polylines, circles and splines as i-arrays. */
//=============================================================================
E_Int K_IO::GenIO::xfigread(
  char* file, char*& varString, E_Int NptsCurve, E_Int NptsLine,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames)
{
  // Discretisation : vectoriel -> arrays
  E_Float pi = 4*atan(1.);
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: xfigread: cannot open file %s.\n", file);
    return 1;
  }

  // Lecture de l'entete
  E_Int res;
  E_Float val;
  E_Float xp0, yp0, zp0, xp1, yp1, zp1;
  char* buf = new char[BUFSIZE];
  res = readGivenKeyword(ptrFile, "#FIG");
  if (res == 0)
  {
    printf("Warning: xfigread: cannot find FIG version.\n");
    fclose(ptrFile);
    return 1;
  }
  res = readDouble(ptrFile, val, -1);
  if (res != 0 && val < 3.2)
  {
    printf("Warning: xfigread: can only read files with version greater than 3.2.\n");
    fclose(ptrFile);
    return 1;
  }

  // paper position (Landscape,...)
  res = readWord(ptrFile, buf); 

  // disposition (Center,...)
  res = readWord(ptrFile, buf); 

  // Unites
  res = readWord(ptrFile, buf); 

  // Paper format (Letter,...)
  res = readWord(ptrFile, buf); 

  // Zoom
  res = readDouble(ptrFile, val, -1);

  // Precision (Single, Double)
  res = readWord(ptrFile, buf); 
  
  // ?? 
  res = readDouble(ptrFile, val, -1);
   
  // Scale
  E_Float scale;
  res = readDouble(ptrFile, scale, -1);
  res = readDouble(ptrFile, val, -1);
  scale = K_FUNC::E_max(scale, 1.e-6);
  //printf("scale " SF_F_ "\n", scale);

  // Lecture des elements
  varString = new char [8];
  strcpy(varString, "x,y,z");
  E_Int type = 2;
  res = readInt(ptrFile, type, -1);

  FldArrayF* an = NULL;
  E_Float depth = 0.;

  while (res != 2)
  {
    // type
    switch (type)
    {
      case 1: // ellipse
      {
        E_Int sousType = 0;
        res = readInt(ptrFile, sousType, -1);
        // sous-type 1 : ellipse centre + rayon
        for (E_Int i = 0; i < 10; i++)
        {
          res = readDouble(ptrFile, val, -1);
          if (i == 4) depth = val;
        }  
        
        res = readDouble(ptrFile, val, -1);
        E_Float xr1 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float yr1 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float a = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float b = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float xr = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float yr = val/scale;
        res = readDouble(ptrFile, val, -1);
        res = readDouble(ptrFile, val, -1);
        an = new FldArrayF(NptsCurve, 3);
        FldArrayF& coord = *an;
        
        if (sousType == 2 || sousType == 4)
        { xr = xr1; yr = yr1;}

        for (E_Int i = 0; i < NptsCurve; i++)
        {
          E_Float t = i*2*pi/(NptsCurve-1);
          coord(i,1) = xr + a*cos(t);
          coord(i,2) = yr + b*sin(t);
          coord(i,3) = depth;
        }

        structField.push_back(an);
        ni.push_back(NptsCurve);
        nj.push_back(1);
        nk.push_back(1);
      }
      break;
        
      case 2: // poly line
      {
        E_Int sousType = 0;
        res = readInt(ptrFile, sousType, -1);
        for (E_Int i = 0; i < 13; i++)
        {
          res = readDouble(ptrFile, val, -1);
          if (i == 4) depth = val;
        }
        E_Int npts = 0;
        res = readInt(ptrFile, npts, -1);

        E_Int nt = (npts-1)*(NptsLine-1)+1;
        an = new FldArrayF(nt, 3);
        FldArrayF& coord = *an;

        res = readDouble(ptrFile, val, -1);
        xp0 = val/scale;
        res = readDouble(ptrFile, val, -1);
        yp0 = val/scale;
        zp0 = depth;
        coord(0,1) = xp0;
        coord(0,2) = yp0;
        coord(0,3) = zp0;
        E_Float alpha = 1. / (NptsLine-1);
        
        E_Int nc = 1;
        for (E_Int i = 1; i < npts; i++)
        {
          res = readDouble(ptrFile, val, -1);
          xp1 = val/scale;
          res = readDouble(ptrFile, val, -1);
          yp1 = val/scale;
          zp1 = depth;
          for (E_Int j = 1; j < NptsLine; j++)
          {
            coord(nc, 1) = xp0 + j*alpha*(xp1-xp0);
            coord(nc, 2) = yp0 + j*alpha*(yp1-yp0);
            coord(nc, 3) = zp0 + j*alpha*(zp1-zp0); nc++;
          }
          xp0 = xp1; yp0 = yp1; zp0 = zp1;
        }
        if (sousType == 1) // polyline ouverte
        {
          structField.push_back(an);
          ni.push_back(nt);
          nj.push_back(1);
          nk.push_back(1);
        }
        else
        {
          FldArrayI* cn = new FldArrayI(nt-1, 2);
          FldArrayI& c = *cn;
          for (E_Int i = 0; i < nt-1; i++)
          {
            c(i,1) = i+1;
            c(i,2) = i+2;
          }
          c(nt-2,2) = 1;
          unstructField.push_back(an);
          connect.push_back(cn);
          eltType.push_back(1); // BAR
        }
      }
      break;

      case 3: // spline
      {
        E_Int sousType = 0;
        res = readInt(ptrFile, sousType, -1);

        for (E_Int i = 0; i < 11; i++)
        {
          res = readDouble(ptrFile, val, -1);
          if (i == 4) depth = val;
        }
        E_Int npts = 0;
        res = readInt(ptrFile, npts, -1);
     
        FldArrayF coord(npts, 3);
 
        // control points
        for (E_Int i = 0; i < npts; i++)
        {
          res = readDouble(ptrFile, val, -1);
          coord(i,1) = val/scale;

          res = readDouble(ptrFile, val, -1);
          coord(i,2) = val/scale;
 
          coord(i,3) = depth;
        }

        // control weights?
        for (E_Int i = 0; i <= npts; i++)
        {
          res = readDouble(ptrFile, val, -1); 
        }

        E_Int order = 3;

        FldArrayF* PF = new FldArrayF(NptsCurve,3);
        K_COMPGEOM::spline(npts, order, NptsCurve,
                           coord.begin(1),coord.begin(2),coord.begin(3),*PF);

        structField.push_back(PF);
        ni.push_back(NptsCurve);
        nj.push_back(1);
        nk.push_back(1);
        
      }
      break;

      case 4: // texte
      {
        printf("Warning: xfigread: text discarded...\n");
        E_Int sousType = 0;
        res = readInt(ptrFile, sousType, -1);
        for (E_Int i = 0; i < 11; i++)
        {
          res = readDouble(ptrFile, val, -1);
        }
        int l = 1;
        while (l > 0)
        {
          res = readWord(ptrFile, buf);
          l = strlen(buf);
          if (l > 3 && buf[l-1] == '1' && buf[l-2] == '0' 
              && buf[l-3] == '0' && buf[l-4] == '\\')
            l = -1;
        }
      }
      break;

      case 5: // portion de cercle
      {
        E_Int sousType = 0;
        res = readInt(ptrFile, sousType, -1);
        for (E_Int i = 0; i < 12; i++)
        {
          res = readDouble(ptrFile, val, -1);
          if (i == 4) depth = val;
        }
        res = readDouble(ptrFile, val, -1);
        E_Float xr = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float yr = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float x1 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float y1 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float x2 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float y2 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float x3 = val/scale;
        res = readDouble(ptrFile, val, -1);
        E_Float y3 = val/scale;

        // On cherche le cercle sous la forme x = xr + cos(at+b), 
        // y = yr + sin(at+b), passant par P1 (t=0) et P3 (t=1)
        //printf("PR " SF_F2_ "\n", xr, yr);
        //printf("P1 " SF_F2_ "\n", x1, y1);
        //printf("P2 " SF_F2_ "\n", x2, y2);
        //printf("P3 " SF_F2_ "\n", x3, y3);
        
        // sens du parcours :
        E_Float prod = (x1-xr)*(y2-yr)-(x2-xr)*(y1-yr);
        //printf("prod = " SF_F_ "\n", prod);
        E_Float a, b, R, t;
        R = (x1-xr)*(x1-xr)+(y1-yr)*(y1-yr);
        R = sqrt(R);
        if (K_FUNC::fEqualZero(x1-xr) == true)
          b = -pi/2.;
        else
          b = atan((y1-yr)/(x1-xr));
        
        if (K_FUNC::fEqualZero(xr + R*cos(b)-x1, 1.e-6) != true || 
            K_FUNC::fEqualZero(yr + R*sin(b)-y1, 1.e-6) != true)
          b = b +  pi;

        if (K_FUNC::fEqualZero(x3-xr) == true)
          a = -pi/2. - b;
        else
          a = atan((y3-yr)/(x3-xr)) - b;

        if (K_FUNC::fEqualZero(xr + R*cos(a+b)-x3, 1.e-6) != true || 
            K_FUNC::fEqualZero(yr + R*sin(a+b)-y3, 1.e-6) != true)
          a = a +  pi;

        if (a > 0 && prod < 0) a = -2*pi + a;
        if (a < 0 && prod > 0) a = 2*pi + a;

        an = new FldArrayF(NptsCurve, 3);
        FldArrayF& coord = *an;
        
        for (E_Int i = 0; i < NptsCurve; i++)
        {
          t = i*1./(NptsCurve-1);
          coord(i,1) = xr + R*cos(a*t+b);
          coord(i,2) = yr + R*sin(a*t+b);
          coord(i,3) = depth;
        }

        structField.push_back(an);
        ni.push_back(NptsCurve);
        nj.push_back(1);
        nk.push_back(1);
      }
      break;

      default: ;
    }
    res = readDouble(ptrFile, val, -1); type = E_Int(val);
  }

  // Cree les noms des zones structurees
  E_Int structSize = structField.size();
  for (E_Int i = 0; i < structSize; i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, i);
    zoneNames.push_back(zoneName);
  }

  // Cree les noms des zones non structurees
  for (unsigned int i=0; i < unstructField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone" SF_D_, i+structSize);
    zoneNames.push_back(zoneName);
  }

  delete [] buf;
  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Write arrays as xfig polyLines.
   Others are discarded. */
//=============================================================================
E_Int K_IO::GenIO::xfigwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  // All zones must have posx, posy, posz
  E_Int posx, posy, posz, ind;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: xfigwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Ecriture de l'entete
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: xfigwrite: I can't open file %s.\n", file);
    return 1;
  }

  E_Float scale = 1200.;
  fprintf(ptrFile, "#FIG 3.2\nLandscape\nCenter\nInches\nLetter\n100.00\nSingle\n-2\n1200 2\n");
  
  // Calcul du depth : on met le z moyen de chaque zone
  E_Int nzones = structField.size();
  E_Int nzoneu = unstructField.size();
  E_Int nzonet = nzones + nzoneu;
  FldArrayI depthi(nzonet);
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    FldArrayF& f = *structField[zone];
    E_Float depth = 0.;
    for (E_Int i = 0; i < f.getSize(); i++)
      depth = depth + f(i, posz);
    depth = depth / f.getSize();
    depthi[zone] = E_Int(depth);
  }
  for (E_Int zone = 0; zone < nzoneu; zone++)
  {
    FldArrayF& f = *unstructField[zone];
    E_Float depth = 0.;
    for (E_Int i = 0; i < f.getSize(); i++)
      depth = depth + f(i, posz);
    depth = depth / f.getSize();
    depthi[zone+nzones] = E_Int(depth);
  }

  // On rescale depth entre 0 et 999
  E_Int depthMin = 1000000;
  E_Int depthMax = -1000000;
  for (E_Int zone = 0; zone < nzonet; zone++)
  {
    depthMin = K_FUNC::E_min(depthMin, depthi[zone]);
    depthMax = K_FUNC::E_max(depthMax, depthi[zone]);
  }
  if (depthMin == depthMax)
    for (E_Int zone = 0; zone < nzonet; zone++)
      depthi[zone] = 50;
  else
  {
    for (E_Int zone = 0; zone < nzonet; zone++)
      depthi[zone] = E_Int((depthi[zone]-depthMin)*999./(depthMax - depthMin));
  }
  
  // Structured -> open polylines for each constant lines
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    E_Int nil = ni[zone];
    E_Int njl = nj[zone];
    E_Int nkl = nk[zone];
    FldArrayF& f = *structField[zone];

    for (E_Int k = 0; k < nkl; k++)
      for (E_Int j = 0; j < njl; j++)
      {
        fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 " SF_D_ "\n", 
                depthi[zone], nil);
        fprintf(ptrFile, "\t");
        for (E_Int i = 0; i < nil; i++)
        {
          ind = i + j*nil + k*nil*njl;
          fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                  E_Int(f(ind,posy)*scale));
        }
        fprintf(ptrFile, "\n");
      }

    if (njl != 1)
    {
      for (E_Int k = 0; k < nkl; k++)
        for (E_Int i = 0; i < nil; i++)
        {
          fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 " SF_D_ "\n", 
                  depthi[zone], njl);
          fprintf(ptrFile, "\t");
          for (E_Int j = 0; j < njl; j++)
          {
            ind = i + j*nil + k*nil*njl;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
          }
          fprintf(ptrFile, "\n");
        }
    }

    if (nkl != 1)
    {
      for (E_Int j = 0; j < njl; j++)
        for (E_Int i = 0; i < nil; i++)
        {
          fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 " SF_D_ "\n", 
                  depthi[zone], nkl);
          fprintf(ptrFile, "\t");
          for (E_Int k = 0; k < nkl; k++)
          {
            ind = i + j*nil + k*nil*njl;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
          }
          fprintf(ptrFile, "\n");
        }
    }
  }

  // Unstructured -> closed polylines
  for (E_Int zone = 0; zone < nzoneu; zone++)
  {
    FldArrayF& f = *unstructField[zone];
    FldArrayI* cm = connect[zone];
    E_Int nc = cm->getNConnect();

    for (E_Int n = 0; n < nc; n++)
    {
      FldArrayI& c = *(cm->getConnect(n));
      E_Int elt = eltTypes[zone][n];   
      if (elt != 1 && elt != 2 && elt != 3 && elt != 4 && elt != 6 && elt != 7)
      {                       
        printf("Error: fig: unrecognised element type: %d.\n", elt);
        continue;
      }

      // for each element
      for (E_Int i = 0; i < c.getSize(); i++)
      {
        switch (elt)
        {
          case 1: // BAR
            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 2);
            fprintf(ptrFile, "\t");
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");
            break;

          case 2: // TRI
            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 4);
            fprintf(ptrFile, "\t");
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,1) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");
            break;
          
          case 3: // QUAD
            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 5);
            fprintf(ptrFile, "\t");
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,4) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,1) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");
            break;

          case 4: // TETRA
            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 4);
            fprintf(ptrFile, "\t");
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 3);
            fprintf(ptrFile, "\t");
            ind = c(i,1) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,4) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 2);
            fprintf(ptrFile, "\t");
            ind = c(i,4) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");
            break;

          case 6: // PENTA (PRISM)
            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 4);
            fprintf(ptrFile, "\t");
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 4);
            fprintf(ptrFile, "\t");
            ind = c(i,1) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,4) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,5) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 3);
            fprintf(ptrFile, "\t");
            ind = c(i,4) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,6) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 2);
            fprintf(ptrFile, "\t");
            ind = c(i,6) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,5) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");
            break;

          case 7: // HEXA
            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 5);
            fprintf(ptrFile, "\t");
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,4) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,1) - 1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 4);
            fprintf(ptrFile, "\t");
            ind = c(i,1) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,5) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,6) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,2) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 4);
            fprintf(ptrFile, "\t");
            ind = c(i,4) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,8) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,7) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,3) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 2);
            fprintf(ptrFile, "\t");
            ind = c(i,5) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,8) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");

            fprintf(ptrFile, "2 1 0 1 0 7 " SF_D_ " -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    depthi[zone+nzones], 2);
            fprintf(ptrFile, "\t");
            ind = c(i,6) -1;
            fprintf(ptrFile, SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            ind = c(i,7) -1;
            fprintf(ptrFile,SF_D2_ " ", E_Int(f(ind,posx)*scale), 
                    E_Int(f(ind,posy)*scale));
            fprintf(ptrFile, "\n");
            break;

          default:
            printf("Error: xfig: unrecognised element type: %d.\n", elt);
        }
      }
    }
  }

  fclose(ptrFile);
  return 0;
}
