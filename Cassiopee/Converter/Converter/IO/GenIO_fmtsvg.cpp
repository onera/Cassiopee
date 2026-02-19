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

// Formated svg file support
  
#define UPDATESTACK \
  rx0 = stackX[0]; ry0 = stackY[0]; \
  for (size_t i = 1; i < stackX.size(); i++) \
  { rx0 += stackX[i]; ry0 += stackY[i]; stackX[i] = rx0; stackY[i] = ry0;}

#define UPDATESTACK2 \
  E_Int ncurves = (npts-1)/3; \
  for (E_Int j = 0; j < ncurves; j++) \
  { rx0 = stackX[3*j]; ry0 = stackY[3*j]; \
  for (E_Int i = 1; i < 4; i++) \
  { stackX[i+3*j] += rx0; stackY[i+3*j] += ry0;}}

#define PRINTSTACK2 ;
#define PRINTSTACK \
  printf("stack:\n"); \
  for (size_t i = 0; i < stackX.size(); i++) printf( SF_D_ ": %g %g\n",i,stackX[i],stackY[i]);
    
#define GENERATE \
  if (commandP == 0 || commandP == 1 || commandP == 2 || commandP == 3) \
    { \
      if (commandP == 1) { UPDATESTACK } \
      else if (commandP == 3) { UPDATESTACK } \
      else if (commandP == 5) { UPDATESTACK } \
      PRINTSTACK2; \
      if (command == 6) an = new FldArrayF(npts+1, 3); \
      else an = new FldArrayF(npts, 3); \
      FldArrayF& coord = *an; \
      for (E_Int i = 0; i < npts; i++) \
      { \
        coord(i, 1) = stackX[i]; \
        coord(i, 2) = -stackY[i]; \
        coord(i, 3) = 0.; \
      } \
      if (command == 6) {coord(npts,1) = coord(0,1); coord(npts,2) = coord(0,2); coord(npts,3) = coord(0,3);} \
      structField.push_back(an); \
      ni.push_back(coord.getSize()); nj.push_back(1); nk.push_back(1); \
    } \
    else if (commandP == 4 || commandP == 5) \
    { \
      if (commandP == 5) { UPDATESTACK2 } \
      PRINTSTACK2; \
      FldArrayF bezierPts(4, 3); \
      E_Int ncurves = (npts-1)/3; \
      for (E_Int j = 0; j < ncurves; j++) \
      { \
        for (E_Int i = 0; i < 4; i++) \
          { \
            bezierPts(i, 1) = stackX[i+j*3]; \
            bezierPts(i, 2) = -stackY[i+j*3]; \
            bezierPts(i, 3) = 0.; \
          } \
          if (density == -1. && exportCtrlPts == false) \
          { \
            an = new FldArrayF(NptsCurve, 3); \
            K_COMPGEOM::bezier(bezierPts.getSize(), NptsCurve, \
                               bezierPts.begin(1), \
                               bezierPts.begin(2), \
                               bezierPts.begin(3), *an); \
          } \
          else if (exportCtrlPts == true) { \
            an = new FldArrayF(4, 3); \
            (*an) = bezierPts; \
          } \
          else { \
            FldArrayF tmp; \
            K_COMPGEOM::regularBezier(bezierPts.getSize(), -1, density, \
                                      bezierPts.begin(1), \
                                      bezierPts.begin(2), \
                                      bezierPts.begin(3), tmp);\
            an = new FldArrayF(tmp.getSize(),3); \
            *an = tmp; \
          } \
          structField.push_back(an); \
          ni.push_back(an->getSize()); nj.push_back(1); nk.push_back(1); \
        }\
    }

# include <stdio.h>
# include <stdlib.h>
# include <vector>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "Def/DefFunction.h"
# include "Connect/connect.h"
# include "CompGeom/compGeom.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* svgread 
   Read svg polyLines and Bezier curves as i-arrays.
   On peut regler le nombre de pts par courbe 
   (si density=-1) ou la densite. */
//=============================================================================
E_Int K_IO::GenIO::svgread(
  char* file, char*& varString, E_Float density,
  E_Int NptsCurve, E_Int NptsLine, 
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames)
{
  bool exportCtrlPts = false;
  if (NptsCurve == 0) { exportCtrlPts = true; NptsCurve=100; }
  FldArrayF bezierPts(4, 3); // pts de controle Bezier cubique
  FldArrayF tmp(NptsCurve, 3);
  
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: svgread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");

  // Lecture de l'entete
  E_Int res;
  char* buf = new char[BUFSIZE];
  char* prevData = new char[BUFSIZE];
  char* keyword = new char[BUFSIZE];
  E_Float zlayer = 0.;
  E_Float stack[3];
  vector<E_Float> stackX;
  vector<E_Float> stackY;
  E_Int command = 0; E_Int commandP = 0;
  E_Float rx0, ry0;

  res = readGivenKeyword(ptrFile, "<?");
  if (res == 0)
  {
    printf("Warning: svgread: cannot find xml header.\n");
    fclose(ptrFile);
    return 1;
  }
  res = readGivenKeyword(ptrFile, ">");
  if (res == 0)
  {
    printf("Warning: svgread: file seems to be corrupted.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Lecture du fichier
  list<const char*> knownKeywords;
  knownKeywords.push_back("SODIPODI:TYPE");
  knownKeywords.push_back("SODIPODI:CX");
  knownKeywords.push_back("SODIPODI:CY");
  knownKeywords.push_back("SODIPODI:R1");
  knownKeywords.push_back("SODIPODI:R2");
  knownKeywords.push_back("SODIPODI:ARG1");
  knownKeywords.push_back("SODIPODI:ARG2");
  knownKeywords.push_back("INKSCAPE:FLATSIDE");
  knownKeywords.push_back("INKSCAPE:ROUNDED");
  knownKeywords.push_back("INKSCAPE:RANDOMIZED");
  knownKeywords.push_back("SODIPODI:EXPANSION");
  knownKeywords.push_back("SODIPODI:REVOLUTION");
  knownKeywords.push_back("SODIPODI:RADIUS");
  knownKeywords.push_back("SODIPODI:ARGUMENT");
  knownKeywords.push_back("SODIPODI:T0");
  knownKeywords.push_back("STYLE");
  knownKeywords.push_back("ID");
  knownKeywords.push_back("D");
  knownKeywords.push_back("TRANSFORM");

  res = readWord(ptrFile, buf);
  while (res >= 0)
  {
    // Commentaire
    if (strcmp(buf, "<!") == 0)
    {
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }
    // Entete svg
    else if (strcmp(buf, "<svg") == 0)
    {
      //printf("Entete svg\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Entete sopodi
    else if (strcmp(buf, "<sopodi:namedview") == 0)
    {
      //printf("Entete sopodi\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Entete metadata
    else if (strcmp(buf, "<metadata") == 0)
    {
      //printf("Bloc metadata\n");
      res = readGivenKeyword(ptrFile, "/METADATA>");
      if (res == 0) res = -1;
    }
    
    // Defs
    else if (strcmp(buf, "<defs") == 0)
    {
      //printf("Bloc defs\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Graphic layer
    else if (strcmp(buf, "<g") == 0)
    {
      zlayer = zlayer + 1.;
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Text
    else if (strcmp(buf, "<text") == 0)
    {
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Path => Objet
    else if (strcmp(buf, "<path") == 0)
    {
      FldArrayF* an = NULL;
      res = readDataAndKeyword(ptrFile, buf, knownKeywords, prevData, keyword);
      if (res == 1) res = -1;
      E_Bool found = false;
      while (res >= 0 && found == false) // OK
      {
        if (strcmp(keyword, "D") == 0)
        {
          res = readWord(ptrFile, buf);
          E_Int reread = true;
          stackX.clear(); stackY.clear();
          command = 0; commandP = -1;
          while (res >= 0)
          {
            reread = true;
            //printf(">> buf=%s\n", buf);
            E_Int ls = strlen(buf);
            E_Int lastChar = buf[ls-1];
            if (lastChar == 'M') command = 0; // absolute move
            else if (lastChar == 'm')  // relative move
            { command = 1; }
            else if (lastChar == 'L') // absolute line 
            { command = 2; } 
            else if (lastChar == 'l') // relative line 
            { command = 3; } 
            else if (lastChar == 'C') // absolute bezier 3 
            { command = 4; } 
            else if (lastChar == 'c') // relative bezier 3 
            { command = 5; } 
            else if (lastChar == 'Z') // close path 
            { command = 6; } 
            else if (lastChar == 'z') // close path 
            { command = 6; } 
            else if (lastChar == 'Q') // absolute bezier 2 
            { command = 7; } 
            else if (lastChar == 'q') // relative bezier 2 
            { command = 8; }
            else if (ls >= 2 && buf[ls-2] == 'z')
            { command = 6; reread=false; }

            if (commandP != command) // nouvelle commande
            {
              if (commandP != -1) // il faut executer commandP
              {
                E_Int npts = stackX.size();
                if (npts > 1)
                {
                  GENERATE;
                }
                // Conserve le dernier pt pour continuer
                stack[0] = stackX.back();
                stack[1] = stackY.back();
                stackX.clear(); stackY.clear();
                stackX.push_back(stack[0]);
                stackY.push_back(stack[1]);
              }
              commandP = command;
            }
            else // not a command, must be a point
            {
              reread = readTwoCoordinates(buf, ptrFile, stack); // 0 = fin de chaine
              stackX.push_back(stack[0]);
              stackY.push_back(stack[1]);
              if (reread == false) // il faut engendrer la courbe suivant command
              {
                E_Int npts = stackX.size();
                if (npts > 1)
                {
                  GENERATE;
                }
                stackX.clear(); stackY.clear();
              }
            }
            if (reread == true) res = readWord(ptrFile, buf);
            else res = -1;
          } // while res >= 0
          found = true;
          res = readGivenKeyword(ptrFile, ">"); // fin du bloc path
          if (res == 0) res = -1;
        }
        else if (strcmp(keyword, "STYLE") == 0)
        {
          res = readGivenKeyword(ptrFile, "\"");
          res = readGivenKeyword(ptrFile, "\"");
          if (res == 0) res = -1;
        }
        else if (strcmp(keyword, "ID") == 0)
        {
          res = readGivenKeyword(ptrFile, "\"");
          res = readGivenKeyword(ptrFile, "\"");
          if (res == 0) res = -1;
        }
        if (found == false)
          res = readDataAndKeyword(ptrFile, buf, 
                                   knownKeywords, prevData, keyword);
        if (res == 1) res = -1;
      } // bloc path
    }
    res = readWord(ptrFile, buf);
  }
  // Cree les noms des zones
  for (size_t i = 0; i < structField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%zu", i);
    zoneNames.push_back(zoneName);
  }

  delete [] buf;
  delete [] prevData;
  delete [] keyword;
  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Write arrays as svg polyLines.
   Others are discarded. */
//=============================================================================
E_Int K_IO::GenIO::svgwrite(
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
    printf("Warning: svgwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Build writing data format
  char format1[82], format2[85], format3[84], format4[85];
  char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format{i}
  sprintf(format1," %s,%s", dataFmtl, dataFmtl);
  sprintf(format2," L %s,%s", dataFmtl, dataFmtl);
  sprintf(format3," %s,%s ", dataFmtl, dataFmtl);
  sprintf(format4," L %s,%s ", dataFmtl, dataFmtl);

  // Ecriture de l'entete
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: svgwrite: I can't open file %s.\n", file);
    return 1;
  }

  fprintf(ptrFile, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
  fprintf(ptrFile, "<!-- Created with Converter (Onera) -->\n");

  // Entete svg
  fprintf(ptrFile, "<svg\n   xmlns:svg=\"http://www.w3.org/2000/svg\"\n   xmlns=\"http://www.w3.org/2000/svg\"\n   version=\"1.0\"\n   width=\"10.\"\n   height=\"10.\"\n   id=\"svg2\">\n");
  fprintf(ptrFile, "<defs\n   id=\"defs4\" />\n");
  fprintf(ptrFile, "<g\n   id=\"layer1\">\n");

  E_Int nzones = structField.size();
  E_Int nzoneu = unstructField.size();
  E_Int nc = 0;

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
        fprintf(ptrFile, "<path\n   d=\"M");
        for (E_Int i = 0; i < nil; i++)
        {
          ind = i + j*nil + k*nil*njl;
          if (i == 0)
            fprintf(ptrFile, format1, f(ind,posx), f(ind,posy));
          else
            fprintf(ptrFile, format2, f(ind,posx), f(ind,posy));
        }
        fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_"\" />\n", nc); nc++;
      }

    if (njl != 1)
    {
      for (E_Int k = 0; k < nkl; k++)
        for (E_Int i = 0; i < nil; i++)
        {
          fprintf(ptrFile, "<path\n   d=\"M");
          for (E_Int j = 0; j < njl; j++)
          {
            ind = i + j*nil + k*nil*njl;
            if (j == 0)
              fprintf(ptrFile,format1, f(ind,posx), f(ind,posy));
            else
              fprintf(ptrFile,format2, f(ind,posx), f(ind,posy));
          }
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
        }
    }

    if (nkl != 1)
    {
      for (E_Int j = 0; j < njl; j++)
        for (E_Int i = 0; i < nil; i++)
        {
          fprintf(ptrFile, "<path\n   d=\"M");
          for (E_Int k = 0; k < nkl; k++)
          {
            ind = i + j*nil + k*nil*njl;
            if (k == 0)
              fprintf(ptrFile,format1, f(ind,posx), f(ind,posy));
            else
              fprintf(ptrFile,format2, f(ind,posx), f(ind,posy));
          }
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
        }
    }
  }

  // Unstructured -> closed polylines
  for (E_Int zone = 0; zone < nzoneu; zone++)
  {
    FldArrayF& f = *unstructField[zone];
    FldArrayI* cm = connect[zone];
    E_Int nco = cm->getNConnect();
    
    for (E_Int n = 0; n < nco; n++)
    {
      FldArrayI& c = *(cm->getConnect(n));
      E_Int elt = eltTypes[zone][n];
  
      if (elt != 1 && elt != 2 && elt != 3 && elt != 4 && elt != 6 && elt != 7)
      {                       
        printf("Error: svg: unrecognised element type: " SF_D_ ".\n", elt);
        continue;
      }

      for (E_Int i = 0; i < c.getSize(); i++)
      {
        switch (elt)
        {
          case 1: // BAR
            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) - 1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,2) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
            break;

          case 2: // TRI
            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) - 1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,2) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,1) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
            break;
          
          case 3: // QUAD
            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) - 1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,2) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,4) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,1) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
            break;

          case 4: // TETRA
            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) - 1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,2) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,1) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,4) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,2) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,4) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,3) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
            break;

          case 6: // PENTA (PRISM)
            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) - 1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,2) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,1) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,4) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,5) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,2) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,4) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,6) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,6) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,5) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
            break;

          case 7: // HEXA
            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) - 1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,2) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,4) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,1) - 1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,1) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,5) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,6) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,2) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,4) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,8) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,7) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            ind = c(i,3) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,5) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,8) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;

            fprintf(ptrFile, "<path\n   d=\"M");
            ind = c(i,6) -1;
            fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
            ind = c(i,7) -1;
            fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
            fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path" SF_D_ "\" />\n", nc); nc++;
            break;

          default:
            printf("Error: svg: unrecognised element type: " SF_D_ ".\n", elt);
        }
      }
    }
  }

  fprintf(ptrFile, "   </g>\n</svg>\n");
  fclose(ptrFile);
  return 0;
}
