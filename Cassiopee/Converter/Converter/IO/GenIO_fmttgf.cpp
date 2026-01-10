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

// Formated TGF (TGrid file - Ansys) file support (not usable)

# include <stdio.h>
# include <stdlib.h>
# include <vector>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include "Def/DefFunction.h"
# include "Connect/connect.h"
# include "CompGeom/compGeom.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* tgfread
   Read zone coordinates.
   Read BC faces.
   Retourne NGONs. */
//=============================================================================
E_Int K_IO::GenIO::tgfread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames,
  vector<FldArrayI*>& BCFaces,
  vector<char*>& BCNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: tgfread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");

  // Recherche le nombre de domaines
  E_Int ret, val, noz, id, start, end, type, left, right, i1, i2, i3, dtype;
  E_Float x, y, z;
  E_Int nzones = 0;

  // Read block
  ret = K_IO::GenIO::readGivenKeyword(ptrFile, "(");
  while (ret == 1)
  {
    ret = K_IO::GenIO::readInt(ptrFile, val);
    printf("reading a block " SF_D_ "\n", val);
    switch (val)
    {
      case 0: // comment
        break; 

      case 1: // title
        break;
        
      case 2: // ??
        break;

      case 4: // ??
      {
        K_IO::GenIO::readGivenKeyword(ptrFile, "(");
        K_IO::GenIO::readGivenKeyword(ptrFile, ")");
      }
      break; 

      case 10: // Nodes
      {
        K_IO::GenIO::readGivenKeyword(ptrFile, "(");
        K_IO::GenIO::readHexaInt(ptrFile, id);
        if (id != 0) nzones++;
        K_IO::GenIO::readHexaInt(ptrFile, start);
        K_IO::GenIO::readHexaInt(ptrFile, end);
        K_IO::GenIO::readHexaInt(ptrFile, type);
        K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        if (id != 0) // read coords
        {
          E_Int nvertex = end-start+1;
          FldArrayF* fp = new FldArrayF(nvertex, 3);
          unstructField.push_back(fp);
          //connect.push_back(NULL);
          //eltType.push_back(0);
          char* zoneName = new char[BUFSIZE+1];
          strcpy(zoneName, "zone");
          zoneNames.push_back(zoneName);
          E_Float* xp = fp->begin(1); 
          E_Float* yp = fp->begin(2);
          E_Float* zp = fp->begin(3);
          K_IO::GenIO::readGivenKeyword(ptrFile, "(");
          for (E_Int i = 0; i < nvertex; i++)
          {
            K_IO::GenIO::readDouble(ptrFile, x);
            K_IO::GenIO::readDouble(ptrFile, y);
            K_IO::GenIO::readDouble(ptrFile, z);
            xp[i] = x; yp[i] = y; zp[i] = z;
          }
          K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        }
      }
      break;
      
      case 11: // edges
      {
        K_IO::GenIO::readGivenKeyword(ptrFile, "(");
        K_IO::GenIO::readHexaInt(ptrFile, id);
        K_IO::GenIO::readHexaInt(ptrFile, start);
        K_IO::GenIO::readHexaInt(ptrFile, end);
        K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        if (id != 0)
        {
          FldArrayI* cn = new FldArrayI(end-start+1, 2);
          E_Int* cn1 = cn->begin(1);
          E_Int* cn2 = cn->begin(2);
          //connect.push_back(cn);
          //eltType.push_back(1); // BAR
          K_IO::GenIO::readGivenKeyword(ptrFile, "(");
          for (E_Int i = 0; i < end-start+1; i++)
          {
            K_IO::GenIO::readHexaInt(ptrFile, i1);
            K_IO::GenIO::readHexaInt(ptrFile, i2);
            cn1[i] = i1; cn2[i] = i2;
          }
          K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        }
      }
      break;

      case 12: // cells
      {
        K_IO::GenIO::readGivenKeyword(ptrFile, "(");
        K_IO::GenIO::readGivenKeyword(ptrFile, ")");
      }
      break; 
      
      case 13: // faces
      {
        K_IO::GenIO::readGivenKeyword(ptrFile, "(");
        K_IO::GenIO::readHexaInt(ptrFile, id);
        K_IO::GenIO::readHexaInt(ptrFile, start);
        K_IO::GenIO::readHexaInt(ptrFile, end);
        K_IO::GenIO::readHexaInt(ptrFile, type);
        K_IO::GenIO::readHexaInt(ptrFile, dtype);
        K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        if (id != 0)
        {
          FldArrayI* cn = new FldArrayI(end-start+1, 3);
          E_Int* cn1 = cn->begin(1);
          E_Int* cn2 = cn->begin(2);
          E_Int* cn3 = cn->begin(3);
          if (connect.size() == 0) connect.push_back(cn);
          eltType.push_back(2); // TRI
          K_IO::GenIO::readGivenKeyword(ptrFile, "(");
          for (E_Int i = 0; i < end-start+1; i++)
          {
            K_IO::GenIO::readHexaInt(ptrFile, i1); // TRI
            K_IO::GenIO::readHexaInt(ptrFile, i2);
            K_IO::GenIO::readHexaInt(ptrFile, i3);
            K_IO::GenIO::readHexaInt(ptrFile, right);
            K_IO::GenIO::readHexaInt(ptrFile, left);
            cn1[i] = i1; cn2[i] = i2; cn3[i] = i3;
          }
          K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        }
      }
      break;

      case 38:
      {

      }
      break;

      default:
      {
        K_IO::GenIO::readGivenKeyword(ptrFile, "(");
        K_IO::GenIO::readHexaInt(ptrFile, noz);
        K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        if (noz != 0)
        {
          K_IO::GenIO::readGivenKeyword(ptrFile, "(");
          K_IO::GenIO::readGivenKeyword(ptrFile, ")");
        }
      }
      break;
    }

    K_IO::GenIO::readGivenKeyword(ptrFile, ")"); // end block
    ret = K_IO::GenIO::readGivenKeyword(ptrFile, "("); // next block
  }
    

 //  // Lecture du Nodes
 //  if (val == 10)
 //  {
 //    // read Nodes header
 //    K_IO::GenIO::readGivenKeyword(ptrFile, "(");
 //    K_IO::GenIO::readInt(ptrFile, noz); 
 //    K_IO::GenIO::readWord(ptrFile, buf);
 //    first =  strtol(number, NULL, 16);
 //    K_IO::GenIO::readInt(ptrFile, buf);
 //    last =  strtol(number, NULL, 16);
 //    K_IO::GenIO::readInt(ptrFile, buf);
 //    type =  strtol(number, NULL, 16); // 0: virtual, 1 (any), 2 (boundary)
 //    K_IO::GenIO::readGivenKeyword(ptrFile, ")");
 //    // read nodes
 //    K_IO::GenIO::readGivenKeyword(ptrFile, "(");
 //    for (E_Int i = first; i < last; i++)
 //    {
 //    K_IO::GenIO::readDouble(ptrFile, x); 
 //    K_IO::GenIO::readDouble(ptrFile, y); 
 //    K_IO::GenIO::readDouble(ptrFile, z);
 //    }
 //    K_IO::GenIO::readGivenKeyword(ptrFile, ")");
 //    K_IO::GenIO::readGivenKeyword(ptrFile, ")");
 //  }

 //  // Lecture Faces
 //  if (val == 13)
 // {
 //    // read Faces header
 //    K_IO::GenIO::readGivenKeyword(ptrFile, "(");
 //    K_IO::GenIO::readInt(ptrFile, noz); 
 //    K_IO::GenIO::readWord(ptrFile, buf);
 //    first =  strtol(number, NULL, 16);
 //    K_IO::GenIO::readInt(ptrFile, buf);
 //    last =  strtol(number, NULL, 16);
 //    K_IO::GenIO::readInt(ptrFile, buf);
 //    type =  strtol(number, NULL, 16);
 //    K_IO::GenIO::readInt(ptrFile, buf);
 //    Elttype =  strtol(number, NULL, 16);
 //    K_IO::GenIO::readGivenKeyword(ptrFile, ")");
 //    // read faces
 //    if (eltType == 2) // BAR
 //    {
 //      for (E_Int i = first; i < last; i++)
 //      {
 //        K_IO::GenIO::readInt(ptrFile, n1); 
 //        K_IO::GenIO::readInt(ptrFile, n2);
 //        K_IO::GenIO::readInt(ptrFile, cr); 
 //        K_IO::GenIO::readInt(ptrFile, cl);
 //      }
 //    }
 // }


 //  // Lecture du nombre de domaines
 //  E_Int ret;
 //  E_Int nbDomains, nvertex, nfaces, ncells, nboundaryfaces;
 //  char buf[BUFSIZE+1];
 //  E_Int value, np;
 //  E_Float x, y, z;
 //  ret = readInt(ptrFile, nbDomains);
 //  if (ret == -1) return 1;
 //  if (ret == 1) skipLine(ptrFile);

 //  for (E_Int nb = 0; nb < nbDomains; nb++)
 //  {
 //    skipLine(ptrFile); // ----
 //    skipLine(ptrFile); // 1.0 DONNEES GENERALES
  
 //    ret = readWord(ptrFile, buf); // nom du domaine
 //    char* zoneName = new char[BUFSIZE+1];
 //    strcpy(zoneName, buf);
 //    zoneNames.push_back(zoneName);
 //    if (ret == 2) skipLine(ptrFile);
    
 //    ret = readWord(ptrFile, buf); // type
 //    if (ret == 2) skipLine(ptrFile);
 //    if (strcmp(buf, "3D") != 0) // 2D ou autre
 //      skipLine(ptrFile); // skip definitions 2D
    
 //    skipLine(ptrFile); // echelle de longueur
 //    ret = readInt(ptrFile, nvertex);
 //    if (ret == 1) skipLine(ptrFile);
 //    //printf("nbre de vertex=" SF_D_ "\n", nvertex);

 //    ret = readInt(ptrFile, nfaces); 
 //    if (ret == 1) skipLine(ptrFile);
 //    //printf("nbre de faces=" SF_D_ "\n", nfaces);

 //    ret = readInt(ptrFile, ncells); 
 //    if (ret == 1) skipLine(ptrFile);
 //    //printf("nbre de cellules=" SF_D_ "\n", ncells);
    
 //    ret = readInt(ptrFile, nboundaryfaces); 
 //    if (ret == 1) skipLine(ptrFile);

 //    FldArrayF* fp = new FldArrayF(nvertex, 3);
 //    unstructField.push_back(fp);
 //    E_Float* xp = fp->begin(1); 
 //    E_Float* yp = fp->begin(2);
 //    E_Float* zp = fp->begin(3);
 //    eltType.push_back(8);
    
 //    // Coordonnees
 //    skipLine(ptrFile);
 //    for (E_Int i = 0; i < nvertex; i++)
 //    {
 //      readInt(ptrFile, value);
 //      readDouble(ptrFile, x); 
 //      readDouble(ptrFile, y); 
 //      ret = readDouble(ptrFile, z);
 //      xp[i] = x; yp[i] = y; zp[i] = z;
 //      //printf(SF_F3_ "\n", x,y,z);
 //    }
 //    if (ret == 1) skipLine(ptrFile);

 //    // Connectivite Faces->Noeuds
 //    // On est oblige de la lire en storage variable
 //    FldArrayI connect1(nfaces*10);
 //    E_Int size = 0;
 //    skipLine(ptrFile);
    
 //    for (E_Int i = 0; i < nfaces; i++)
 //    {
 //      readInt(ptrFile, value); // no de la face
 //      ret = readInt(ptrFile, np); connect1[size] = np; size++;
 //      //printf("np=" SF_D_ "\n", np);
 //      if (size+np+1 > connect1.getSize()) 
 //        connect1.reAlloc(size+(nfaces-i+10)*(np+1));
 //      for (E_Int j = 0; j < np; j++)
 //      {
 //        ret = readInt(ptrFile, value); connect1[size] = value; size++;
 //      }
 //    }
 //    connect1.reAlloc(size);
 //    //printf(SF_D_ "\n", size);

 //    // Connectivite faces->elts
 //    skipLine(ptrFile);
 //    FldArrayI connect2(nfaces*2);
 //    E_Int* cn2 = connect2.begin();
 //    for (E_Int i = 0; i < nfaces; i++)
 //    {
 //      readInt(ptrFile, value); //printf("no=" SF_D_ "\n", value); // no face
 //      readInt(ptrFile, value);  //printf("nv=" SF_D_ "\n", value); // nbre de faces valides (1 ou 2)
 //      ret = readInt(ptrFile, value); cn2[2*i] = value;
 //      ret = readInt(ptrFile, value); cn2[2*i+1] = value;
 //      //printf(SF_D2_ "\n", cn2[2*i], cn2[2*i+1]);
 //    }
 //    if (ret == 1) skipLine(ptrFile);

 //    // Construit la connectivite elts->faces
 //    // Compte les faces pour chaque elements
 //    FldArrayI nbface(ncells); nbface.setAllValuesAtNull();
 //    E_Int* nbfacep = nbface.begin();
 //    E_Int a, b;
 //    for (E_Int i = 0; i < nfaces; i++)
 //    {
 //      a = cn2[2*i]; b = cn2[2*i+1];
 //      nbfacep[a-1]++;
 //      if (b != 0) nbfacep[b-1]++;
 //    }
 //    //for (E_Int i = 0; i < ncells; i++) printf(SF_D2_ "\n", i, nbfacep[i]);

 //    size = ncells;
 //    for (E_Int i = 0; i < ncells; i++) size += nbfacep[i];
 //    FldArrayI pos(ncells); // position des faces pour chaque elt
 //    E_Int* posp = pos.begin();
 //    posp[0] = 0;
 //    for (E_Int i = 0; i < ncells-1; i++)
 //    {
 //      posp[i+1] = posp[i]+(nbfacep[i]+1);
 //    }

 //    FldArrayI connect3(size); E_Int* connect3p = connect3.begin();
 //    for (E_Int i = 0; i < ncells; i++)
 //    {
 //      connect3p[posp[i]] = nbface[i];
 //    }

 //    nbface.setAllValuesAtNull();
 //    for (E_Int i = 0; i < nfaces; i++)
 //    {
 //      a = cn2[2*i]; b = cn2[2*i+1];
 //      connect3p[posp[a-1]+nbfacep[a-1]+1] = i+1; nbfacep[a-1]++;
 //      if (b != 0) {connect3p[posp[b-1]+nbfacep[b-1]+1] = i+1; nbfacep[b-1]++;}
 //    }
 //    pos.malloc(0); connect2.malloc(0); nbface.malloc(0);

 //    // Fusionne les connectivites
 //    FldArrayI* cn = new FldArrayI(connect1.getSize()+connect3.getSize()+4, 1);
 //    connect.push_back(cn);
 //    E_Int* cp = cn->begin();
 //    cp[0] = nfaces;
 //    cp[1] = connect1.getSize();    
 //    for (E_Int i = 0; i < connect1.getSize(); i++)
 //    {
 //      cp[i+2] = connect1[i];
 //    }
 //    int pt = 2+connect1.getSize();
 //    cp[pt] = ncells;
 //    cp[pt+1] = connect3.getSize();
 //    for (E_Int i = 0; i < connect3.getSize(); i++)
 //    {
 //      cp[pt+2+i] = connect3[i];
 //    }
 //    connect3.malloc(0); connect1.malloc(0);

 //    // Lit les frontieres
 //    skipLine(ptrFile);
 //    FldArrayI* faces = new FldArrayI(nboundaryfaces);
 //    char* names = new char [nboundaryfaces*BUFSIZE];
 //    E_Int* facesp = faces->begin();
 //    E_Int le;
 //    E_Int c = 0;

 //    for (E_Int i = 0; i < nboundaryfaces; i++)
 //    {
 //      ret = readInt(ptrFile, value); // numerotation (skip) 
 //      ret = readInt(ptrFile, value); // no de la face
 //      facesp[i] = value;
 //      ret = readWord(ptrFile, buf); // nom de la BC
 //      le = strlen(buf);
 //      for (E_Int l = 0; l < le; l++) { names[c+l] = buf[l]; }
 //      c += le;
 //      names[c] = '\0'; c++;
 //    }
    
 //    BCFaces.push_back(faces);
 //    BCNames.push_back(names);
    
 //  } // Blocs

  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Write arrays as tgf file.
   Write only NGON_n arrays. */
//=============================================================================
E_Int K_IO::GenIO::tgfwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames,
  PyObject* BCFaces)
{
  printf("tgfwrite: not implemented.\n");
  return 1;
}
