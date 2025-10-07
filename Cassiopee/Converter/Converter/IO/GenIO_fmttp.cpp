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
// Formated tecplot file support

# include <stdio.h>
# include <string.h>

#include <stdlib.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "../converter.h"

using namespace K_ARRAY;
using namespace std;
using namespace K_FLD;

//const E_Int K_IO::BUFSIZE = 8192;
//const E_Int K_IO::BUFSIZE = 512;

//=============================================================================
/* tpread */
//=============================================================================
E_Int K_IO::GenIO::tpread(
  char* file, char*& varString,
  std::vector<FldArrayF*>& structField,
  std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
  std::vector<FldArrayF*>& unstructField,
  std::vector<FldArrayI*>& connectivity,
  std::vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: tpread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [K_ARRAY::VARSTRINGLENGTH];

  char buf[BUFSIZE+1];
  char keyword[BUFSIZE+1];
  char nextKeyword[BUFSIZE+1];
  char prevData[BUFSIZE+1];
  char data[BUFSIZE+1];
  list<const char*> knownKeywords;
  E_Int c, l, i, j;
  E_Int nvar = 0;
  E_Int sizet = 0;
  E_Int packing = 0;
  E_Int zone = 0;
  E_Int nil, njl, nkl;
  E_Int np = 0; E_Int ne = 0; E_Int nfaces = 0;
  E_Int sizeElt = 0;
  E_LONG pos = 0;
  E_Int nheadlines=0;
  char lastSeparator='\0';
  /* Local vector for structured and unstructured zones names */
  vector<char*> structZoneNames, unstructZoneNames;
  char zoneName[BUFSIZE+1];
  E_Int ZONETFound, VARIABLESFound;

  knownKeywords.push_back("TITLE");
  knownKeywords.push_back("TEXT");
  knownKeywords.push_back("VARIABLES");
  knownKeywords.push_back("VARLOCATION"); // dummy
  knownKeywords.push_back("STRANDID");
  knownKeywords.push_back("SOLUTIONTIME");
  knownKeywords.push_back("DATAPACKING");
  knownKeywords.push_back("PARENTZONE");
  knownKeywords.push_back("ZONETYPE");
  knownKeywords.push_back("FILETYPE"); // dummy
  knownKeywords.push_back("ZONE"); // zone + title
  knownKeywords.push_back("AUXDATA");
  knownKeywords.push_back("FACENODES");
  knownKeywords.push_back("NODES");
  knownKeywords.push_back("ELEMENTS"); // nbre d'elements
  knownKeywords.push_back("BOUNDARYFACES");
  knownKeywords.push_back("BOUNDARYCONNECTIONS");
  knownKeywords.push_back("FACES"); // nbre de faces
  knownKeywords.push_back("DT");
  knownKeywords.push_back("ET"); // type d'element
  knownKeywords.push_back("CS"); // dummy Text Position
  knownKeywords.push_back("HU"); // dummy Text Position
  knownKeywords.push_back(" I"); // ni
  knownKeywords.push_back(",I"); // ni
  knownKeywords.push_back(" J"); // nj
  knownKeywords.push_back(",J"); // nj
  knownKeywords.push_back(" K"); // nk
  knownKeywords.push_back(",K"); // nk
  knownKeywords.push_back(" F"); // format (point, block)
  knownKeywords.push_back(",F"); // format (point, block)
  knownKeywords.push_back(" N"); // nbre de pts
  knownKeywords.push_back(",N"); // nbre de pts
  knownKeywords.push_back(" E"); // elements
  knownKeywords.push_back(",E"); // elements
  knownKeywords.push_back(" Y"); // dummy
  knownKeywords.push_back(",Y"); // dummy
  knownKeywords.push_back(" T"); // dummy Zone Title
  knownKeywords.push_back(",T"); // dummy Zone Title
  knownKeywords.push_back(" X"); // dummy Text Position
  knownKeywords.push_back(",X"); // dummy Text Position
  knownKeywords.push_back(" H"); // dummy Text Position
  knownKeywords.push_back(",H"); // dummy Text Position

  /* Header read */
  E_Int res = readKeyword(ptrFile, keyword);

  blocread: ;
  np = 0; ne = 0; nil = 0; njl = 0; nkl = 0; nfaces = 0;
  // must be unique, otherwise short zone
  ZONETFound = 0; VARIABLESFound = 0;
  strcpy(data, keyword);

  while (res == 0) // boucle sur les keywords
  {
    compressString(keyword);

    nheadlines++;
    pos = KFTELL(ptrFile);
    res = readDataAndKeyword(ptrFile, buf,
                             knownKeywords,
                             prevData, nextKeyword);

    // Pour une zone avec des petites data, le ZONE_T, VARIABLES ou
    // DATAPACKING a deja ete trouve, mais les data non lues
    if (strcmp(nextKeyword, "TEXT") == 0) res=1;
    if (strcmp(nextKeyword, "ZONE") == 0 && ZONETFound == 1) res = 1;
    if (strcmp(nextKeyword, "VARIABLES") == 0 && VARIABLESFound == 1) res = 1;

    if (res == 1)
    {
      // on tronque prevData avant les data float
      l = BUFSIZE;
      for (i = 0; i < l; i++)
      {
        if (prevData[i] == '\n')
        {
          prevData[i] = '\0'; break;
        }
      }
      for (j = i+1; j < l; j++) data[j-i-1] = prevData[j];
      data[l-i-1] = '\0';
      strcpy(nextKeyword, "");
    }

    if (strcmp(keyword, "VARIABLES") == 0)
    {
      VARIABLESFound = 1;
      // Build varString from prevData
      //printf("Variables: %s", prevData);
      l = strlen(prevData);
      nvar = 0;
      c = 0;
      for (i = 0; i < l; i++)
      {
        if (prevData[i] == '"' || prevData[i] == '\'')
        {
          if (lastSeparator=='\0')
          {
            lastSeparator=prevData[i];
            nvar++;
          }
          else 
          {
            lastSeparator='\0';
            prevData[i]=' ';
            if (c != 0) {varString[c] = ','; c++;}
          }
        }
        else if ((prevData[i] == ' ' || prevData[i] == '\n' || prevData[i] == ','
            || prevData[i] == '\r') && lastSeparator=='\0')
        {
          if (i > 0 && prevData[i-1] == ' ') ;
          else if (i > 0 && prevData[i-1] == '\n') ;
          else if (i > 0 && prevData[i-1] == '\r') ;
          else if (i > 0 && prevData[i-1] == ',') ;
          else if (i == 0) ;
          else
          {
            nvar++;
            if (c != 0) {varString[c] = ','; c++;}
          }
        }
        else
        {
          if (prevData[i] != '"')
          {
            varString[c] = prevData[i]; c++;
          }
        }
      }
      if (c > 0) varString[c-1] = '\0'; // final ,
      else varString[c] = '\0';
    }
    else if (strcmp(keyword, "ZONE") == 0)
    {
      ZONETFound = 1;
      zoneName[0] = '\0';
    }
    else if (strcmp(keyword, "T") == 0 || strcmp(keyword, ",T") == 0)
    {
      // Build zoneName from prevData
      c = 0;
      l = strlen(prevData);
      for (i = 0; i < l; i++)
      {
        if ((prevData[i] != ' ')&&(prevData[i] != '"')&&
            (prevData[i] != ',')&&(prevData[i] != '\n')&& c<BUFSIZE)
        {
          zoneName[c] = prevData[i]; c++;
        }
      }
      zoneName[c] = '\0';
      //printf("Trouver: %s\n", zoneName);
    }
    else if (strcmp(keyword, "ET") == 0 || strcmp(keyword, "ZONETYPE") == 0)
    {
      if (matchInString(prevData, "LINESEG") == 1)
      { eltType.push_back(1); sizeElt = 2; } // BAR
      else if (matchInString(prevData, "TRIANGLE") == 1)
      { eltType.push_back(2); sizeElt = 3; } // TRIANGLE
      else if (matchInString(prevData, "QUADRILATERAL") == 1)
      { eltType.push_back(3); sizeElt = 4; } // QUAD
      else if (matchInString(prevData, "TETRAHEDRON") == 1)
      { eltType.push_back(4); sizeElt = 4; } // TETRA
      else if (matchInString(prevData, "BRICK") == 1)
      { eltType.push_back(7); sizeElt = 8; } // HEXA
      else if (matchInString(prevData, "FEPOLYHEDRON") == 1)
      { eltType.push_back(8); sizeElt = -1; } // NGON volumique
      else if (matchInString(prevData, "FEPOLYGON") == 1)
      { eltType.push_back(8); sizeElt = -2; } // NGON surfacique
      else if (matchInString(prevData, "ORDERED") == 1)
      { ; }
      else
      {printf("Warning: tpread: unknown element type.\n");
        fclose(ptrFile); return 1;}
    }
    else if (strcmp(keyword, "FACES") == 0)
    {
      nfaces = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "I") == 0 || strcmp(keyword, ",I") == 0)
    {
      //printf("i: %s\n", prevData);
      nil = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "J") == 0 || strcmp(keyword, ",J") == 0)
    {
      //printf("j: %s\n", prevData);
      njl = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "K") == 0 || strcmp(keyword, ",K") == 0)
    {
      //printf("k: %s\n", prevData);
      nkl = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "N") == 0 || strcmp(keyword, "NODES") == 0 || strcmp(keyword, ",N") == 0)
    {
      np = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "E") == 0 || strcmp(keyword, "ELEMENTS") == 0 || strcmp(keyword, ",E") == 0)
    {
      ne = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "F") == 0 || strcmp(keyword, ",F") == 0)
    {
      compressString(prevData);
      //printf("format : %s\n", prevData);
      if (matchInString(prevData, "BLOCK") == 1 ||
          matchInString(prevData, "FEBLOCK") == 1)
        packing = 0;
      else packing = 1;
    }
    else if (strcmp(keyword, "DATAPACKING") == 0)
    {
      compressString(prevData);
      if (matchInString(prevData, "BLOCK") == 1) packing = 0;
      else packing = 1;
    }
    strcpy(keyword, nextKeyword);
    for (size_t i = 0; i < strlen(prevData); i++)
      if (prevData[i] == '\n') nheadlines++;
  }

  // Allocating arrays
  if (nil == 0 && np == 0)
  {
    // Try to guess nil from file content
    KFSEEK(ptrFile, pos, SEEK_SET);
    // On compte le nombre de ligne de donnees du fichier
    E_Int ndatalines = -1; // On part de la dernière ligne de header (Faux si aucun header)
    char firstline[BUFSIZE+1];
    do {
      i = 0;
      buf[i] = '\0';
      do { // read one line
        c = fgetc(ptrFile);
        if (c != EOF) buf[i++] = (char)c;
      } while(c != EOF && c != '\n');
      buf[i]='\0';
      if (ndatalines == 0) 
      {
        strcpy(firstline,buf);
      } 
      if ('\n' == c) {
        ++ndatalines;
      } else if (c==EOF) { // Verifie que la derniere ligne est vide sinon ajoute une ligne
        i=0;
        while(buf[i]!='\0') {
          if (buf[i++]!=' ') {
            ndatalines++;
            break;
          }
        }
      }
    } while(EOF != c);
    // On compte le nombre de colonne de la premiere ligne
    i = 0;
    char ch;
    E_Int ncol=1;
    do 
    { // read one line
      do 
      {
        ch = firstline[i++];
      } while(ch == ' ');
      if (firstline[i] == ' ') ncol++;
    } while(ch != '\n' && ch != '\0');
    if (nvar == 0) 
    { // Essaie de lire un fichier sans entete (gnuplot/tecplot fmt point)
      nvar = ncol; // Format point par defaut
      ndatalines++;
      for (i=0; i < nvar; i++)
      {
        varString[i*3]='V';
        varString[i*3+1]='1'+i; //static_cast<char>(i);
        varString[i*3+2]=',';
      }
      varString[strlen(varString)-1]='\0';
    }
    if (ncol == nvar) 
    { // Format point
      nil = ndatalines;
      packing = 1;
    } 
    else 
    {
      if (ndatalines==nvar) 
      { // Format block simple
        nil=ncol;
        packing=0;
      }
    }
    if (nil == 0) 
    {
      printf("Warning: tpread: can not read ni,nj,nk or np.\n");
      fclose(ptrFile);
      return 1;
    }
  }
  if (nvar == 0)
  {
    printf("Warning: tpread: can not read variables.\n");
    fclose(ptrFile);
    return 1;
  }

  if (nil != 0)
  { ni.push_back(nil);
    if (njl != 0) nj.push_back(njl);
    else nj.push_back(1);
    if (nkl != 0) nk.push_back(nkl);
    else nk.push_back(1);
  }

  if (nil != 0) sizet = ni[zone]*nj[zone]*nk[zone];
  else sizet = np;

  FldArrayF* f = new FldArrayF(sizet, nvar);
  if (nil != 0) structField.push_back(f);
  else unstructField.push_back(f);
  if (nil != 0)
  {
    // put zone name in structured list
    E_Int structSize = structZoneNames.size();
    if (zoneName[0] == '\0') 
    sprintf(zoneName, "StructZone" SF_D_, structSize);
    char* name = new char[BUFSIZE+1]; strcpy(name, zoneName);
    structZoneNames.push_back(name);
    zoneName[0] = '\0';
  }
  else
  {
    // put zone name in unstructured list
    E_Int unstructSize = unstructZoneNames.size();
    if (zoneName[0] == '\0')
    sprintf(zoneName, "UnstructZone" SF_D_, unstructSize);
    char* name = new char[BUFSIZE+1]; strcpy(name, zoneName);
    unstructZoneNames.push_back(name);
    zoneName[0] = '\0';
  }

  // Reading data: data contains the begining of read data
  KFSEEK(ptrFile, pos, SEEK_SET); skipLine(ptrFile);

  E_Float value;
  E_Int ret = 0;
  c = 0;
  E_Float* ff = f->begin();
  E_Int size = sizet*nvar;

  if (packing == 0)
  {
    // Lecture format block
    for (i = 0; i < size; i++)
    {
      ret = readDouble(ptrFile, value);
      if (ret == -1)
      {
        printf("Warning: tpread: can not read a block of data.\n");
        fclose(ptrFile);
        return 1;
      }
      ff[i] = value;
    }
  }
  else
  {
    // Lecture format point
    for (i = 0; i < sizet; i++)
    {
      for (j = 0; j < nvar; j++)
      {
        ret = readDouble(ptrFile, value);
        if (ret == -1)
        {
          printf("Warning: tpread: can not read enough points.\n");
          fclose(ptrFile);
          return 1;
        }
        ff[i+j*sizet] = value;
      }
    }
  }

  // Lecture connectivite (eventuellement)
  if (np != 0)
  {
    E_Int index;
    if (sizeElt >= 0) // Basic elements
    {
      FldArrayI* cn = new FldArrayI(ne, sizeElt);
      connectivity.push_back(cn);
      for (E_Int i = 0; i < ne; i++)
      {
        for (E_Int n = 1; n <= sizeElt; n++)
        {
          ret = readInt(ptrFile, index);
          (*cn)(i, n) = index;
        }
      }
    }
    else // NGON
    {
      // Read node count par face (info locale)
      FldArrayI count(nfaces);
      E_Int* countp = count.begin();
      E_Int size = 0;
      if (sizeElt == -1) // NGon volumique uniquement
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
          ret = readInt(ptrFile, index); countp[i] = index;
          size += index+1;
        }
      }
      else // surfacique
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
           countp[i] = 2; // toujours 2 en surfacique
           size += 3;
        }
      }

      // Read face nodes
      FldArrayI faces(size);
      E_Int* ptr = faces.begin();
      E_Int n;

      for (E_Int i = 0; i < nfaces; i++)
      {
        n = countp[i];
        ptr[0] = n;
        for (E_Int j = 0; j < n; j++)
        {
          ret = readInt(ptrFile, index);
          ptr[j+1] = index;
        }
        ptr += n+1;
      }

      // connectivite FE
      FldArrayI cFE(nfaces, 2);
      E_Int* cFE1 = cFE.begin(1); E_Int* cFE2 = cFE.begin(2);
      for (E_Int i = 0; i < nfaces; i++)
      {
        ret = readInt(ptrFile, index);
        cFE1[i] = index;
      }
      for (E_Int i = 0; i < nfaces; i++)
      {
        ret = readInt(ptrFile, index);
        cFE2[i] = index;
      }

      // Construit la connectivite elts->faces
      FldArrayI cEF;
      K_CONNECT::connectFE2EF(cFE, ne, cEF);
      // Fusionne les connectivites
      FldArrayI* cn = new FldArrayI(faces.getSize()+cEF.getSize()+4);
      connectivity.push_back(cn);
      E_Int* cp = cn->begin();
      cp[0] = nfaces;
      cp[1] = faces.getSize();
      for (E_Int i = 0; i < faces.getSize(); i++)
      {
        cp[i+2] = faces[i];
      }
      int pt = 2+faces.getSize();
      cp[pt] = ne;
      cp[pt+1] = cEF.getSize();
      for (E_Int i = 0; i < cEF.getSize(); i++)
      {
        cp[pt+2+i] = cEF[i];
      }
    }
  }
  // }
  while (1)
  {
    // next bloc
    pos = KFTELL(ptrFile);
    res = readKeyword(ptrFile, keyword);

    if (res!=0) break;
    else
    {
      compressString(keyword);
      if (strcmp(keyword, "ZONE")==0) break;
      else
      {
        KFSEEK(ptrFile, pos, SEEK_SET);
        skipLine(ptrFile);
      }
    }
  }

  if (strcmp(keyword, "ZONE")==0) {zone++; KFSEEK(ptrFile, pos, SEEK_SET); goto blocread;}

  // Concatenation of structured and unstructured zones names lists
  zoneNames = structZoneNames;
  zoneNames.insert(zoneNames.end(), unstructZoneNames.begin(),
                   unstructZoneNames.end());
  fclose(ptrFile);
  return 0;
}

E_Int K_IO::GenIO::tpread(
  char* file, char*& varString,
  std::vector<FldArrayF*>& structField,
  std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
  std::vector<FldArrayF*>& unstructField,
  std::vector<FldArrayI*>& connectivity,
  std::vector<std::vector<E_Int> >& eltType, vector<char*>& zoneNames,
  E_Int api)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: tpread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [K_ARRAY::VARSTRINGLENGTH];
  E_Int neltTypes = 0;
  eltType.resize(1);

  char buf[BUFSIZE+1];
  char keyword[BUFSIZE+1];
  char nextKeyword[BUFSIZE+1];
  char prevData[BUFSIZE+1];
  char data[BUFSIZE+1];
  list<const char*> knownKeywords;
  E_Int c, l, i, j;
  E_Int nvar = 0;
  E_Int sizet = 0;
  E_Int packing = 0;
  size_t zone = 0, uzone = 0;
  E_Int nil, njl, nkl;
  E_Int np = 0; E_Int ne = 0; E_Int nfaces = 0;
  E_Int ngonDim = -1;
  E_LONG pos = 0;
  E_Int nheadlines=0;
  char lastSeparator='\0';
  /* Local vector for structured and unstructured zones names */
  vector<char*> structZoneNames, unstructZoneNames;
  char zoneName[BUFSIZE+1];
  E_Int ZONETFound, VARIABLESFound;

  knownKeywords.push_back("TITLE");
  knownKeywords.push_back("TEXT");
  knownKeywords.push_back("VARIABLES");
  knownKeywords.push_back("VARLOCATION"); // dummy
  knownKeywords.push_back("STRANDID");
  knownKeywords.push_back("SOLUTIONTIME");
  knownKeywords.push_back("DATAPACKING");
  knownKeywords.push_back("PARENTZONE");
  knownKeywords.push_back("ZONETYPE");
  knownKeywords.push_back("FILETYPE"); // dummy
  knownKeywords.push_back("ZONE"); // zone + title
  knownKeywords.push_back("AUXDATA");
  knownKeywords.push_back("FACENODES");
  knownKeywords.push_back("NODES");
  knownKeywords.push_back("ELEMENTS"); // nbre d'elements
  knownKeywords.push_back("BOUNDARYFACES");
  knownKeywords.push_back("BOUNDARYCONNECTIONS");
  knownKeywords.push_back("FACES"); // nbre de faces
  knownKeywords.push_back("DT");
  knownKeywords.push_back("ET"); // type d'element
  knownKeywords.push_back("CS"); // dummy Text Position
  knownKeywords.push_back("HU"); // dummy Text Position
  knownKeywords.push_back(" I"); // ni
  knownKeywords.push_back(",I"); // ni
  knownKeywords.push_back(" J"); // nj
  knownKeywords.push_back(",J"); // nj
  knownKeywords.push_back(" K"); // nk
  knownKeywords.push_back(",K"); // nk
  knownKeywords.push_back(" F"); // format (point, block)
  knownKeywords.push_back(",F"); // format (point, block)
  knownKeywords.push_back(" N"); // nbre de pts
  knownKeywords.push_back(",N"); // nbre de pts
  knownKeywords.push_back(" E"); // elements
  knownKeywords.push_back(",E"); // elements
  knownKeywords.push_back(" Y"); // dummy
  knownKeywords.push_back(",Y"); // dummy
  knownKeywords.push_back(" T"); // dummy Zone Title
  knownKeywords.push_back(",T"); // dummy Zone Title
  knownKeywords.push_back(" X"); // dummy Text Position
  knownKeywords.push_back(",X"); // dummy Text Position
  knownKeywords.push_back(" H"); // dummy Text Position
  knownKeywords.push_back(",H"); // dummy Text Position

  /* Header read */
  E_Int res = readKeyword(ptrFile, keyword);

  blocread: ;
  np = 0; ne = 0; nil = 0; njl = 0; nkl = 0; nfaces = 0;
  // must be unique, otherwise short zone
  ZONETFound = 0; VARIABLESFound = 0;
  strcpy(data, keyword);

  while (res == 0) // boucle sur les keywords
  {
    compressString(keyword);

    nheadlines++;
    pos = KFTELL(ptrFile);
    res = readDataAndKeyword(ptrFile, buf,
                             knownKeywords,
                             prevData, nextKeyword);

    // Pour une zone avec des petites data, le ZONE_T, VARIABLES ou
    // DATAPACKING a deja ete trouve, mais les data non lues
    if (strcmp(nextKeyword, "TEXT") == 0) res=1;
    if (strcmp(nextKeyword, "ZONE") == 0 && ZONETFound == 1) res = 1;
    if (strcmp(nextKeyword, "VARIABLES") == 0 && VARIABLESFound == 1) res = 1;

    if (res == 1)
    {
      // on tronque prevData avant les data float
      l = BUFSIZE;
      for (i = 0; i < l; i++)
      {
        if (prevData[i] == '\n')
        {
          prevData[i] = '\0'; break;
        }
      }
      for (j = i+1; j < l; j++) data[j-i-1] = prevData[j];
      data[l-i-1] = '\0';
      strcpy(nextKeyword, "");
    }

    if (strcmp(keyword, "VARIABLES") == 0)
    {
      VARIABLESFound = 1;
      // Build varString from prevData
      //printf("Variables: %s", prevData);
      l = strlen(prevData);
      nvar = 0;
      c = 0;
      for (i = 0; i < l; i++)
      {
        if (prevData[i] == '"' || prevData[i] == '\'')
        {
          if (lastSeparator=='\0')
          {
            lastSeparator=prevData[i];
            nvar++;
          }
          else 
          {
            lastSeparator='\0';
            prevData[i]=' ';
            if (c != 0) {varString[c] = ','; c++;}
          }
        }
        else if ((prevData[i] == ' ' || prevData[i] == '\n' || prevData[i] == ','
               || prevData[i] == '\r') && lastSeparator=='\0')
        {
          if (i > 0 && prevData[i-1] == ' ') ;
          else if (i > 0 && prevData[i-1] == '\n') ;
          else if (i > 0 && prevData[i-1] == '\r') ;
          else if (i > 0 && prevData[i-1] == ',') ;
          else if (i == 0) ;
          else
          {
            nvar++;
            if (c != 0) {varString[c] = ','; c++;}
          }
        }
        else
        {
          if (prevData[i] != '"')
          {
            varString[c] = prevData[i]; c++;
          }
        }
      }
      if (c > 0) varString[c-1] = '\0'; // final ,
      else varString[c] = '\0';
    }
    else if (strcmp(keyword, "ZONE") == 0)
    {
      ZONETFound = 1;
      zoneName[0] = '\0';
    }
    else if (strcmp(keyword, "T") == 0 || strcmp(keyword, ",T") == 0)
    {
      // Build zoneName from prevData
      c = 0;
      l = strlen(prevData);
      for (i = 0; i < l; i++)
      {
        if ((prevData[i] != ' ')&&(prevData[i] != '"')&&
            (prevData[i] != ',')&&(prevData[i] != '\n')&& c<BUFSIZE)
        {
          zoneName[c] = prevData[i]; c++;
        }
      }
      zoneName[c] = '\0';
    }
    else if (strcmp(keyword, "ET") == 0 || strcmp(keyword, "ZONETYPE") == 0)
    {
      // In api 1 for basic elements, add the first BE only
      if (matchInString(prevData, "LINESEG") == 1)
      {
        if (api == 3 || neltTypes == 0) eltType[uzone].push_back(1);
        neltTypes++;
      }
      else if (matchInString(prevData, "TRIANGLE") == 1)
      {
        if (api == 3 || neltTypes == 0) eltType[uzone].push_back(2);
        neltTypes++;
      }
      else if (matchInString(prevData, "QUADRILATERAL") == 1)
      {
        if (api == 3 || neltTypes == 0) eltType[uzone].push_back(3);
        neltTypes++;
      }
      else if (matchInString(prevData, "TETRAHEDRON") == 1)
      {
        if (api == 3 || neltTypes == 0) eltType[uzone].push_back(4);
        neltTypes++;
      }
      else if (matchInString(prevData, "BRICK") == 1)
      {
        if (api == 3 || neltTypes == 0) eltType[uzone].push_back(7);
        neltTypes++;
      }
      else if (matchInString(prevData, "FEPOLYHEDRON") == 1) // NGON volumique
      {
        eltType[uzone].push_back(8);
        ngonDim = 3; neltTypes++;
      }
      else if (matchInString(prevData, "FEPOLYGON") == 1) // NGON surfacique
      {
        eltType[uzone].push_back(8);
        ngonDim = 2; neltTypes++;
      }
      else if (matchInString(prevData, "ORDERED") == 1)
      { ; }
      else
      {
        printf("Warning: tpread: unknown element type.\n");
        fclose(ptrFile); return 1;
      }
    }
    else if (strcmp(keyword, "FACES") == 0)
    {
      nfaces = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "I") == 0 || strcmp(keyword, ",I") == 0)
    {
      nil = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "J") == 0 || strcmp(keyword, ",J") == 0)
    {
      njl = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "K") == 0 || strcmp(keyword, ",K") == 0)
    {
      nkl = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "N") == 0 || strcmp(keyword, "NODES") == 0 || strcmp(keyword, ",N") == 0)
    {
      np = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "E") == 0 || strcmp(keyword, "ELEMENTS") == 0 || strcmp(keyword, ",E") == 0)
    {
      ne = convertString2Int(prevData);
    }
    else if (strcmp(keyword, "F") == 0 || strcmp(keyword, ",F") == 0)
    {
      compressString(prevData);
      //printf("format : %s\n", prevData);
      if (matchInString(prevData, "BLOCK") == 1 ||
          matchInString(prevData, "FEBLOCK") == 1)
        packing = 0;
      else packing = 1;
    }
    else if (strcmp(keyword, "DATAPACKING") == 0)
    {
      compressString(prevData);
      if (matchInString(prevData, "BLOCK") == 1) packing = 0;
      else packing = 1;
    }
    strcpy(keyword, nextKeyword);
    for (size_t i = 0; i < strlen(prevData); i++)
      if (prevData[i] == '\n') nheadlines++;
  }

  // Allocating arrays
  if (nil == 0 && np == 0)
  {
    // Try to guess nil from file content
    KFSEEK(ptrFile, pos, SEEK_SET);
    // On compte le nombre de ligne de donnees du fichier
    E_Int ndatalines = -1; // On part de la dernière ligne de header (Faux si aucun header)
    char firstline[BUFSIZE+1];
    do {
      i = 0;
      buf[i] = '\0';
      do { // read one line
        c = fgetc(ptrFile);
        if (c != EOF) buf[i++] = (char)c;
      } while(c != EOF && c != '\n');
      buf[i]='\0';
      if (ndatalines == 0) 
      {
        strcpy(firstline,buf);
      } 
      if ('\n' == c) {
        ++ndatalines;
      } else if (c==EOF) { // Verifie que la derniere ligne est vide sinon ajoute une ligne
        i=0;
        while(buf[i]!='\0') {
          if (buf[i++]!=' ') {
            ndatalines++;
            break;
          }
        }
      }
    } while(EOF != c);
    // On compte le nombre de colonne de la premiere ligne
    i = 0;
    char ch;
    E_Int ncol=1;
    do 
    { // read one line
      do 
      {
        ch = firstline[i++];
      } while(ch == ' ');
      if (firstline[i] == ' ') ncol++;
    } while(ch != '\n' && ch != '\0');
    if (nvar == 0) 
    { // Essaie de lire un fichier sans entete (gnuplot/tecplot fmt point)
      nvar = ncol; // Format point par defaut
      ndatalines++;
      for (i=0; i < nvar; i++)
      {
        varString[i*3]='V';
        varString[i*3+1]='1'+i; //static_cast<char>(i);
        varString[i*3+2]=',';
      }
      varString[strlen(varString)-1]='\0';
    }
    if (ncol == nvar) 
    { // Format point
      nil = ndatalines;
      packing = 1;
    } 
    else 
    {
      if (ndatalines==nvar) 
      { // Format block simple
        nil=ncol;
        packing=0;
      }
    }
    if (nil == 0) 
    {
      printf("Warning: tpread: can not read ni,nj,nk or np.\n");
      fclose(ptrFile);
      return 1;
    }
  }
  if (nvar == 0)
  {
    printf("Warning: tpread: can not read variables.\n");
    fclose(ptrFile);
    return 1;
  }

  if (nil != 0)
  { ni.push_back(nil);
    if (njl != 0) nj.push_back(njl);
    else nj.push_back(1);
    if (nkl != 0) nk.push_back(nkl);
    else nk.push_back(1);
  }

  if (nil != 0) sizet = ni[zone]*nj[zone]*nk[zone];
  else sizet = np;

  FldArrayF* f = new FldArrayF(sizet, nvar);
  if (nil != 0) structField.push_back(f);
  else unstructField.push_back(f);
  if (nil != 0)
  {
    // put zone name in structured list
    E_Int structSize = structZoneNames.size();
    if (zoneName[0] == '\0') 
    sprintf(zoneName, "StructZone" SF_D_, structSize);
    char* name = new char[BUFSIZE+1]; strcpy(name, zoneName);
    structZoneNames.push_back(name);
    zoneName[0] = '\0';
  }
  else
  {
    // put zone name in unstructured list
    E_Int unstructSize = unstructZoneNames.size();
    if (zoneName[0] == '\0')
    sprintf(zoneName, "UnstructZone" SF_D_, unstructSize);
    char* name = new char[BUFSIZE+1]; strcpy(name, zoneName);
    unstructZoneNames.push_back(name);
    zoneName[0] = '\0';
  }

  // Reading data: data contains the beginning of read data
  KFSEEK(ptrFile, pos, SEEK_SET); skipLine(ptrFile);

  E_Float value;
  E_Int ret = 0;
  c = 0;
  E_Float* ff = f->begin();
  E_Int size = sizet*nvar;

  if (packing == 0)
  {
    // Lecture format block
    for (i = 0; i < size; i++)
    {
      ret = readDouble(ptrFile, value);
      if (ret == -1)
      {
        printf("Warning: tpread: can not read a block of data.\n");
        fclose(ptrFile);
        return 1;
      }
      ff[i] = value;
    }
  }
  else
  {
    // Lecture format point
    for (i = 0; i < sizet; i++)
    {
      for (j = 0; j < nvar; j++)
      {
        ret = readDouble(ptrFile, value);
        if (ret == -1)
        {
          printf("Warning: tpread: can not read enough points.\n");
          fclose(ptrFile);
          return 1;
        }
        ff[i+j*sizet] = value;
      }
    }
  }

  // Lecture connectivite (eventuellement)
  if (np != 0)
  {
    E_Int ind, elt;
    if (ngonDim == -1) // Basic elements
    {
      vector<E_Int> nvpe(8);
      nvpe[1] = 2; nvpe[2] = 3; nvpe[3] = 4;
      nvpe[4] = 4; nvpe[7] = 8;

      // Create BE (no ME)
      // NB: BE may have be degenerated elements and be an ME but can only
      // be read as a BE here. Clean connectivity in user script if necessary
      elt = eltType[uzone][0];
      E_Int nvpeElt = nvpe[elt];
      FldArrayI* cn2 = new FldArrayI(ne, nvpeElt);
      for (E_Int i = 0; i < ne; i++)
      {
        for (E_Int n = 1; n <= nvpeElt; n++)
        {
          ret = readInt(ptrFile, ind);
          (*cn2)(i, n) = ind;
        }
      }
      connectivity.push_back(cn2);
    }
    else // NGON
    {
      // Reading in api 1
      // Read node count par face (info locale)
      FldArrayI count(nfaces);
      E_Int* countp = count.begin();
      E_Int size = 0;
      if (ngonDim == 3) // NGon volumique uniquement
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
          ret = readInt(ptrFile, ind); countp[i] = ind;
          size += ind+1;
        }
      }
      else // surfacique
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
           countp[i] = 2; // toujours 2 en surfacique
           size += 2+1;
        }
      }

      // Read face nodes
      FldArrayI faces(size);
      E_Int* ptr = faces.begin();
      E_Int n;

      for (E_Int i = 0; i < nfaces; i++)
      {
        n = countp[i];
        ptr[0] = n;
        for (E_Int j = 0; j < n; j++)
        {
          ret = readInt(ptrFile, ind);
          ptr[j+1] = ind;
        }
        ptr += n+1;
      }

      // connectivite FE
      FldArrayI cFE(nfaces, 2);
      E_Int* cFE1 = cFE.begin(1); E_Int* cFE2 = cFE.begin(2);
      for (E_Int i = 0; i < nfaces; i++)
      {
        ret = readInt(ptrFile, ind);
        cFE1[i] = ind;
      }
      for (E_Int i = 0; i < nfaces; i++)
      {
        ret = readInt(ptrFile, ind);
        cFE2[i] = ind;
      }

      // Construit la connectivite elts->faces
      FldArrayI cEF;
      K_CONNECT::connectFE2EF(cFE, ne, cEF);
      // Fusionne les connectivites
      E_Int ngonType = api; // TODO
      E_Int shift = 1; if (ngonType == 3) shift = 0;
      E_Int sizeFN = faces.getSize(); // size in api 1
      E_Int sizeEF = cEF.getSize(); // size in api 1
      if (ngonType == 3) { sizeFN -= nfaces; sizeEF -= ne; }
      PyObject* tpl = K_ARRAY::buildArray3(3, varString, np, ne, nfaces, 
                                           "NGON", sizeFN, sizeEF, ngonType,
                                           false, api);
      FldArrayF* f2; FldArrayI* cn2;
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      E_Int* ngon2 = cn2->getNGon();
      E_Int* nface2 = cn2->getNFace();
      E_Int *indPG2 = NULL, *indPH2 = NULL; 
      if (ngonType == 2 || ngonType == 3) // set offsets
      {
        indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
        indPG2[0] = 0; indPH2[0] = 0;
      }

      E_Int nv, nf, ind1, ind2;
      E_Int c1 = 0, c2 = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        nv = faces[c2]; ngon2[c1] = nv;
        if (ngonType == 3) indPG2[i+1] = indPG2[i] + nv;
        if ((ngonType == 2) and (i+1 < nfaces)) indPG2[i+1] = indPG2[i] + nv;
        for (E_Int j = 0; j < nv; j++)
        {
          ind1 = c1+j+shift; ind2 = c2+j+1;
          ngon2[ind1] = faces[ind2];
        }
        c1 += nv+shift; c2 += nv+1;
      }

      c1 = 0; c2 = 0;
      for (E_Int i = 0; i < ne; i++)
      {
        nf = cEF[c2]; nface2[c1] = nf;
        if (ngonType == 3) indPH2[i+1] = indPH2[i] + nf;
        if ((ngonType == 2) and (i+1 < ne)) indPH2[i+1] = indPH2[i] + nf;
        for (E_Int j = 0; j < nf; j++)
        {
          ind1 = c1+j+shift; ind2 = c2+j+1;
          nface2[ind1] = cEF[ind2];
        }
        c1 += nf+shift; c2 += nf+1;
      }

      connectivity.push_back(cn2);
      delete f2;
    }
  }

  while (1)
  {
    // next bloc
    pos = KFTELL(ptrFile);
    res = readKeyword(ptrFile, keyword);

    if (res != 0) break;
    else
    {
      compressString(keyword);
      if (strcmp(keyword, "ZONE") == 0) break;
      else
      {
        KFSEEK(ptrFile, pos, SEEK_SET);
        skipLine(ptrFile);
      }
    }
  }

  if (strcmp(keyword, "ZONE") == 0)
  {
    zone++;
    if (unstructZoneNames.size() == eltType.size())
    {
      // Last zone was an unstructured zone
      uzone++; eltType.resize(uzone+1);
      neltTypes = 0;
    }
    KFSEEK(ptrFile, pos, SEEK_SET); goto blocread;
  }

  // Last zone was a structured zone
  if (unstructZoneNames.size() != eltType.size()) eltType.resize(eltType.size()-1);

  // Concatenation of structured and unstructured zones names lists
  zoneNames = structZoneNames;
  zoneNames.insert(zoneNames.end(), unstructZoneNames.begin(),
                   unstructZoneNames.end());
  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* tpwrite */
//=============================================================================
E_Int K_IO::GenIO::tpwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  char t[3024];
  E_Int fieldSize, structSize;
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "w");

  if (ptrFile == NULL)
  {
    printf("Warning: tpwrite: cannot open file %s.\n", file);
    return 1;
  }

  // Build writing data format
  char format1[30], format2[60], format3[90], format4[120], format5[150], format6[180];

  char dataFmtl[29];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt);
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format{i}, i=1,6
  sprintf(format1,"%s\n", dataFmtl);
  sprintf(format2,"%s%s\n", dataFmt, dataFmtl);
  sprintf(format3,"%s%s%s\n", dataFmt, dataFmt, dataFmtl);
  sprintf(format4,"%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmtl);
  sprintf(format5,"%s%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmt,
          dataFmtl);
  sprintf(format6,"%s%s%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmt,
          dataFmt, dataFmtl);

  // Write header
  fprintf(ptrFile, "TITLE = \"Generated by Cassiopee\"\n");
  vector<char*> vars;
  extractVars(varString, vars);
  strcpy(t, "VARIABLES = "); fprintf(ptrFile, "%s", t);
  E_Int varsSize = vars.size();
  for (E_Int i = 0; i < varsSize; i++)
  {
    strcpy(t, "\"");
    strcat(t, vars[i]);
    strcat(t, "\" ");
    fprintf(ptrFile, "%s", t);
    delete [] vars[i];
  }
  strcpy(t, "\n"); fprintf(ptrFile, "%s", t);

  // Write zone (structured)
  fieldSize = structField.size();
  structSize = fieldSize;
  for (E_Int cnt = 0; cnt < fieldSize; cnt++)
  {
    FldArrayF& f = *structField[cnt];
    E_Int nijk = ni[cnt]*nj[cnt]*nk[cnt];
    fprintf(ptrFile, "ZONE T=\"%s\",  I=" SF_D_ ",  J=" SF_D_ ",  K=" SF_D_", F=BLOCK\n",
            zoneNames[cnt], ni[cnt], nj[cnt], nk[cnt]);

    for (E_Int n = 1; n <= f.getNfld(); n++)
    {
      E_Float* fn = f.begin(n);
      E_Int i = 0;
      while (i < nijk)
      {
        if (i + 5 < nijk)
        {
          fprintf(ptrFile, format6,
                  fn[i], fn[i+1], fn[i+2],
                  fn[i+3], fn[i+4], fn[i+5]);
          i += 6;
        }
        else break;
      }
      if (nijk-i == 5)
        fprintf(ptrFile, format5,
                fn[i], fn[i+1], fn[i+2],
                fn[i+3], fn[i+4]);
      else if (nijk - i == 4)
        fprintf(ptrFile, format4,
                fn[i], fn[i+1], fn[i+2], fn[i+3]);
      else if (nijk-i == 3)
        fprintf(ptrFile, format3,
                fn[i], fn[i+1], fn[i+2]);
      else if (nijk-i == 2)
        fprintf(ptrFile, format2, fn[i], fn[i+1]);
      else if (nijk-i == 1)
        fprintf(ptrFile, format1, fn[i]);
    }
  }

  // Write zone (unstructured)
  fieldSize = unstructField.size();
  for (E_Int cnt = 0; cnt < fieldSize; cnt++)
  {
    FldArrayF& f = *unstructField[cnt];
    FldArrayI& c = *connect[cnt];
    E_Int nodes = f.getSize();
    E_Int nfld = f.getNfld();
    E_Int elts = c.getSize();
    switch (eltType[cnt])
    {
      case 1: // BAR
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=LINESEG, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 2: // TRI
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=TRIANGLE, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 3: // QUAD
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=QUADRILATERAL, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 4: // TETRA
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=TETRAHEDRON, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 5: // PYRA - FIX as HEXA
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=BRICK, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 6: // PENTA - FIX as HEXA
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=BRICK, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 7: // HEXA
        fprintf(ptrFile,
                "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=BRICK, F=FEBLOCK\n",
                zoneNames[cnt+structSize], nodes, elts);
        break;
      case 8: // NGON
        {
          E_Int nfaces = c.getNFaces();
          E_Int sizeFN = c.getSizeNGon();
          elts = c.getNElts();
          E_Int* ngon = c.getNGon(); 
          E_Int* indPG = c.getIndPG();
          E_Int nf; c.getFace(0, nf, ngon, indPG);
          if (nf > 2) // volumique
            fprintf(ptrFile,
              "ZONE T=\"%s\", Nodes=" SF_D_ ", Elements=" SF_D_ ", Faces=" SF_D_ ", ZONETYPE=FEPOLYHEDRON\nDATAPACKING=BLOCK\nTotalNumFaceNodes=" SF_D_ ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0\n",
              zoneNames[cnt+structSize], nodes, elts, nfaces, sizeFN-nfaces);
          else
            fprintf(ptrFile,
              "ZONE T=\"%s\", Nodes=" SF_D_ ", Elements=" SF_D_ ", Faces=" SF_D_ ", ZONETYPE=FEPOLYGON\nDATAPACKING=BLOCK\nTotalNumFaceNodes=" SF_D_ ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0\n",
              zoneNames[cnt+structSize], nodes, elts, nfaces, sizeFN-nfaces);
        }
        break;
      default:
        printf("tpwrite: unknown type of element. Skipping zone...\n");
        nfld = 0; elts = 0;
    }

    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Int i = 0;
      E_Float* fn = f.begin(n);
      while (i < nodes)
      {
        if (i + 5 < nodes)
        {
          fprintf(ptrFile, format6,
                  fn[i], fn[i+1], fn[i+2],
                  fn[i+3], fn[i+4], fn[i+5]);
          i = i+6;
        }
        else
          break;
      }
      if (nodes-i == 5)
        fprintf(ptrFile, format5,
                fn[i], fn[i+1], fn[i+2],
                fn[i+3], fn[i+4]);
      else if (nodes-i == 4)
        fprintf(ptrFile, format4,
                fn[i], fn[i+1], fn[i+2], fn[i+3]);
      else if (nodes-i == 3)
        fprintf(ptrFile, format3,
                fn[i], fn[i+1], fn[i+2]);
      else if (nodes-i == 2)
        fprintf(ptrFile, format2, fn[i], fn[i+1]);
      else if (nodes-i == 1)
        fprintf(ptrFile, format1, fn[i]);
    }

    // Write connectivity
    if (eltType[cnt] == 5) // PYRA fix as HEXA
    {
      for (E_Int i = 0; i < elts; i++)
      {
        for (E_Int n = 1; n <= c.getNfld(); n++)
        {
          fprintf(ptrFile, SF_D_ " ", c(i,n));
        }
        fprintf(ptrFile, SF_D3_, c(i,5), c(i,5), c(i,5));
        fprintf(ptrFile, "\n");
      }
    }
    else if (eltType[cnt] == 6) // PENTA fix as HEXA
    {
      for (E_Int i = 0; i < elts; i++)
      {
        fprintf(ptrFile, SF_D8_,
                c(i,1), c(i,2), c(i,2), c(i,3),
                c(i,4), c(i,5), c(i,5), c(i,6));
        fprintf(ptrFile, "\n");
      }

    }
    else if (eltType[cnt] == 8) // NGONS
    {
      E_Int nfaces = c.getNFaces();
      E_Int* ngon = c.getNGon();
      E_Int* indPG = c.getIndPG();
      E_Int col = 0;
      E_Int n; c.getFace(0, n, ngon, indPG);

      // node count per face (seult pour les NGON volumiques)
      if (n > 2)
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
          c.getFace(i, n, ngon, indPG);
          fprintf(ptrFile, " " SF_D_, n); col++;
          if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
        }
        fprintf(ptrFile, "\n");
      }
      // face nodes
      col = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        E_Int* face = c.getFace(i, n, ngon, indPG);
        for (E_Int j = 0; j < n; j++) fprintf(ptrFile, " " SF_D_, face[j]);
        col += n;
        if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
      }
      fprintf(ptrFile, "\n");

      // left elements for each face
      FldArrayI cFE;
      K_CONNECT::connectNG2FE(c, cFE);
      col = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        n = cFE(i,1);
        fprintf(ptrFile, " " SF_D_, n); col++;
        if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
      }
      fprintf(ptrFile, "\n");

      // right elements for each face
      col = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        n = cFE(i,2);
        fprintf(ptrFile, " " SF_D_, n); col++;
        if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
      }
      fprintf(ptrFile, "\n");
    }
    else // Other basic elements
    {
      for (E_Int i = 0; i < elts; i++)
      {
        for (E_Int n = 1; n <= c.getNfld(); n++)
        {
          fprintf(ptrFile, SF_D_ " ", c(i,n));
        }
        fprintf(ptrFile, "\n");
      }
    }
  }

  fclose(ptrFile);
  return 0;
}

E_Int K_IO::GenIO::tpwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<vector<E_Int> >& eltType,
  vector<char*>& zoneNames)
{
  char t[3024];
  E_Int fieldSize, structSize;
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "w");

  if (ptrFile == NULL)
  {
    printf("Warning: tpwrite: cannot open file %s.\n", file);
    return 1;
  }

  // Build writing data format
  char format1[30], format2[60], format3[90], format4[120], format5[150], format6[180];
  char format1s[30];

  char dataFmtl[29];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt);
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format{i}, i=1,6
  sprintf(format1s,"%s ", dataFmtl);
  sprintf(format1,"%s\n", dataFmtl);
  sprintf(format2,"%s%s\n", dataFmt, dataFmtl);
  sprintf(format3,"%s%s%s\n", dataFmt, dataFmt, dataFmtl);
  sprintf(format4,"%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmtl);
  sprintf(format5,"%s%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmt,
          dataFmtl);
  sprintf(format6,"%s%s%s%s%s%s\n", dataFmt, dataFmt, dataFmt, dataFmt,
          dataFmt, dataFmtl);

  // Write header
  fprintf(ptrFile, "TITLE = \"Generated by Cassiopee\"\n");
  vector<char*> vars;
  extractVars(varString, vars);
  strcpy(t, "VARIABLES = "); fprintf(ptrFile, "%s", t);
  E_Int varsSize = vars.size();
  for (E_Int i = 0; i < varsSize; i++)
  {
    strcpy(t, "\"");
    strcat(t, vars[i]);
    strcat(t, "\" ");
    fprintf(ptrFile, "%s", t);
    delete [] vars[i];
  }
  strcpy(t, "\n"); fprintf(ptrFile, "%s", t);

  // Write zone (structured)
  fieldSize = structField.size();
  structSize = fieldSize;
  for (E_Int cnt = 0; cnt < fieldSize; cnt++)
  {
    FldArrayF& f = *structField[cnt];
    E_Int nijk = ni[cnt]*nj[cnt]*nk[cnt];
    fprintf(ptrFile, "ZONE T=\"%s\",  I=" SF_D_ ",  J=" SF_D_ ",  K=" SF_D_", F=BLOCK\n",
            zoneNames[cnt], ni[cnt], nj[cnt], nk[cnt]);

    for (E_Int n = 1; n <= f.getNfld(); n++)
    {
      E_Float* fn = f.begin(n);
      E_Int i = 0;
      while (i < nijk)
      {
        if (i + 5 < nijk)
        {
          fprintf(ptrFile, format6,
                  fn[i], fn[i+1], fn[i+2],
                  fn[i+3], fn[i+4], fn[i+5]);
          i += 6;
        }
        else break;
      }
      if (nijk-i == 5)
        fprintf(ptrFile, format5,
                fn[i], fn[i+1], fn[i+2],
                fn[i+3], fn[i+4]);
      else if (nijk - i == 4)
        fprintf(ptrFile, format4,
                fn[i], fn[i+1], fn[i+2], fn[i+3]);
      else if (nijk-i == 3)
        fprintf(ptrFile, format3,
                fn[i], fn[i+1], fn[i+2]);
      else if (nijk-i == 2)
        fprintf(ptrFile, format2, fn[i], fn[i+1]);
      else if (nijk-i == 1)
        fprintf(ptrFile, format1, fn[i]);
    }
  }

  // Write zone (unstructured)
  const E_Int ncmax = 8;
  E_Int nelts;
  E_Int nunsZones = unstructField.size();

  std::vector<const char*> beetfmt; beetfmt.reserve(ncmax);
  beetfmt.push_back("\0");
  beetfmt.push_back("LINESEG\0");
  beetfmt.push_back("TRIANGLE\0");
  beetfmt.push_back("QUADRILATERAL\0");
  beetfmt.push_back("TETRAHEDRON\0");
  for (int i = 0; i < 3; i++) beetfmt.push_back("BRICK\0");

  for (E_Int cnt = 0; cnt < nunsZones; cnt++)
  {
    FldArrayF& f = *unstructField[cnt];
    FldArrayI& c = *connect[cnt];
    const vector<E_Int>& eltTypeZone = eltType[cnt];
    E_Bool isNGon = eltTypeZone[0] == 8;
    E_Int npts = f.getSize();
    E_Int nfld = f.getNfld();
    E_Int packing = 0; // 0: block, 1: point
    // E_Int dim = 2;
    
    // Write zone header
    if (isNGon) // NGON
    {
      nelts = c.getNElts();
      E_Int nfaces = c.getNFaces();
      E_Int sizeFN = c.getSizeNGon();
      E_Int* ngon = c.getNGon(); 
      E_Int* indPG = c.getIndPG();
      E_Int nf; c.getFace(0, nf, ngon, indPG);
      if (nf > 2) // volumique
        fprintf(ptrFile,
          "ZONE T=\"%s\", Nodes=" SF_D_ ", Elements=" SF_D_ ", Faces=" SF_D_ ", ZONETYPE=FEPOLYHEDRON\nDATAPACKING=BLOCK\nTotalNumFaceNodes=" SF_D_ ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0\n",
          zoneNames[cnt+structSize], npts, nelts, nfaces, sizeFN-nfaces);
      else
        fprintf(ptrFile,
          "ZONE T=\"%s\", Nodes=" SF_D_ ", Elements=" SF_D_ ", Faces=" SF_D_ ", ZONETYPE=FEPOLYGON\nDATAPACKING=BLOCK\nTotalNumFaceNodes=" SF_D_ ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0\n",
          zoneNames[cnt+structSize], npts, nelts, nfaces, sizeFN-nfaces);
    }
    else // BE/ME
    {
      E_Int enc = eltTypeZone.size();
      E_Int nc = c.getNConnect();
      if (nc != enc)
      {
        printf("Error: tpwrite: connectivity size mismatch between connect and "
              "eltType, " SF_D_ " and " SF_D_ ".\n", nc, enc);
        return 1;
      }

      // "FEBLOCK": element connectivity in block format, 1 elttype per zone
      // "FEPOINT": element connectivity can contain several elttypes (any order)
      char ffmt[] = "FEBLOCK";
      char etfmt[] = "QUADRILATERAL"; // used if F = "FEPOINT" only, covers all 2D cases
      for (E_Int cnt2 = 0; cnt2 < nunsZones; cnt2++)
      {
        if (eltTypeZone[0] > 3)
        {
          // dim = 3;
          strcpy(etfmt, "BRICK"); // covers all 3D cases
        }
        if (eltTypeZone.size() > 1)
        {
          strcpy(ffmt, "FEPOINT"); packing = 1; break;
        }
      }

      if (nc == 1) 
      {
        strcpy(etfmt, beetfmt[eltTypeZone[0]]);
        nelts = c.getConnect(0)->getSize();
      }
      else
      {
        nelts = 0;
        for (E_Int ic = 0; ic < nc; ic++)
        {
          FldArrayI& cm = *(c.getConnect(ic));
          nelts += cm.getSize();
        }
      }

      fprintf(ptrFile,
              "ZONE T=\"%s\", N=" SF_D_ ", E=" SF_D_ ", ET=%s, F=%s\n",
              zoneNames[cnt+structSize], npts, nelts, etfmt, ffmt);
    }
    
    // Write coordinates
    if (packing == 0) // by block
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Int i = 0;
        E_Float* fn = f.begin(n);
        while (i < npts)
        {
          if (i + 5 < npts)
          {
            fprintf(ptrFile, format6,
                    fn[i], fn[i+1], fn[i+2],
                    fn[i+3], fn[i+4], fn[i+5]);
            i = i+6;
          }
          else
            break;
        }
        if (npts-i == 5)
          fprintf(ptrFile, format5,
                  fn[i], fn[i+1], fn[i+2],
                  fn[i+3], fn[i+4]);
        else if (npts-i == 4)
          fprintf(ptrFile, format4,
                  fn[i], fn[i+1], fn[i+2], fn[i+3]);
        else if (npts-i == 3)
          fprintf(ptrFile, format3,
                  fn[i], fn[i+1], fn[i+2]);
        else if (npts-i == 2)
          fprintf(ptrFile, format2, fn[i], fn[i+1]);
        else if (npts-i == 1)
          fprintf(ptrFile, format1, fn[i]);
      }
    }
    else // by point
    {
      for (E_Int i = 0; i < npts; i++)
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fn = f.begin(n);
          fprintf(ptrFile, format1s, fn[i]);
        }
        fprintf(ptrFile, "\n");
      }
    }
    
    // Write connectivity
    if (isNGon) // NGON
    {
      E_Int* ngon = c.getNGon();
      E_Int* indPG = c.getIndPG();
      E_Int nfaces = c.getNFaces();
      E_Int n; c.getFace(0, n, ngon, indPG);
      E_Int col = 0;

      // node count per face (seult pour les NGON volumiques)
      if (n > 2)
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
          c.getFace(i, n, ngon, indPG);
          fprintf(ptrFile, " " SF_D_, n); col++;
          if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
        }
        fprintf(ptrFile, "\n");
      }

      // face vertices
      col = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        E_Int* face = c.getFace(i, n, ngon, indPG);
        for (E_Int j = 0; j < n; j++) fprintf(ptrFile, " " SF_D_, face[j]);
        col += n;
        if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
      }
      fprintf(ptrFile, "\n");

      // left elements for each face
      FldArrayI cFE;
      K_CONNECT::connectNG2FE(c, cFE);
      col = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        n = cFE(i,1);
        fprintf(ptrFile, " " SF_D_, n); col++;
        if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
      }
      fprintf(ptrFile, "\n");

      // right elements for each face
      col = 0;
      for (E_Int i = 0; i < nfaces; i++)
      {
        n = cFE(i,2);
        fprintf(ptrFile, " " SF_D_, n); col++;
        if (col > 10) { fprintf(ptrFile, "\n"); col = 0; }
      }
      fprintf(ptrFile, "\n");
    }
    else // BE/ME
    {
      E_Int nc = c.getNConnect();
      // E_Int nmax = 8;
      // if (dim == 3) nmax = 8;
      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(c.getConnect(ic));
        E_Int nelts = cm.getSize();
        E_Int nvpe = cm.getNfld();

        // for (E_Int i = 0; i < nelts; i++)
        // {
        //   for (E_Int n = 1; n <= nvpe; n++)
        //     fprintf(ptrFile, SF_D_ " ", cm(i,n));
        //   if (nc > 1 || eltTypeZone[ic] == 5 || eltTypeZone[ic] == 6)
        //   {
        //     // Add placeholders: integer > ntotpts (all zones)
        //     for (E_Int n = nvpe+1; n <= nmax; n++)
        //       fprintf(ptrFile, SF_D_ " ", cm(i,nvpe)); 
        //   }
        //   fprintf(ptrFile, "\n");
        // }

        if (eltTypeZone[ic] == 5) // PYRA fix as HEXA
        {
          for (E_Int i = 0; i < nelts; i++)
          {
            for (E_Int n = 1; n <= nvpe; n++)
              fprintf(ptrFile, SF_D_ " ", cm(i,n));
            fprintf(ptrFile, SF_D3_ "\n", cm(i,5), cm(i,5), cm(i,5));
          }
        }
        else if (eltTypeZone[ic] == 6) // PENTA fix as HEXA
        {
          for (E_Int i = 0; i < nelts; i++)
            fprintf(ptrFile, SF_D8_ "\n",
                    cm(i,1), cm(i,2), cm(i,2), cm(i,3),
                    cm(i,4), cm(i,5), cm(i,5), cm(i,6));
        }
        else
        {
          if (nc == 1 || eltTypeZone[ic] == 1 || eltTypeZone[ic] == 3 || eltTypeZone[ic] == 7)
          {
            for (E_Int i = 0; i < nelts; i++)
            {
              for (E_Int n = 1; n <= nvpe; n++)
                fprintf(ptrFile, SF_D_ " ", cm(i,n)); // standard
              fprintf(ptrFile, "\n");
            }
          }
          else if (eltTypeZone[ic] == 2) // TRI fix as QUAD
          {
            for (E_Int i = 0; i < nelts; i++)
            {
              for (E_Int n = 1; n <= nvpe; n++)
                fprintf(ptrFile, SF_D_ " ", cm(i,n));
              fprintf(ptrFile, SF_D_ "\n", cm(i,3));
            }
          }
          else if (eltTypeZone[ic] == 4) // TETRA fix as HEXA
          {          
            for (E_Int i = 0; i < nelts; i++)
            {
              fprintf(ptrFile, SF_D4_ " ", cm(i,1), cm(i,2), cm(i,3), cm(i,3));
              fprintf(ptrFile, SF_D4_ "\n", cm(i,4), cm(i,4), cm(i,4), cm(i,4));
            }
          }
        }
      }
    }
  }

  fclose(ptrFile);
  return 0;
}