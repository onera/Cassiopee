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

// Formated cedre (Onera) file support

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
/* cedreread
   Read zone coordinates.
   Read BC faces.
   Retourne NGONs. */
//=============================================================================
E_Int K_IO::GenIO::cedreread(
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
    printf("Warning: cedreread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");

  // Lecture du nombre de domaines
  E_Int ret;
  E_Int nbDomains, nvertex, nfaces, ncells, nboundaryfaces;
  char buf[BUFSIZE+1];
  E_Int value, np;
  E_Float x, y, z;
  ret = readInt(ptrFile, nbDomains);
  if (ret == -1) return 1;
  if (ret == 1) skipLine(ptrFile);

  for (E_Int nb = 0; nb < nbDomains; nb++)
  {
    skipLine(ptrFile); // ----
    skipLine(ptrFile); // 1.0 DONNEES GENERALES
  
    ret = readWord(ptrFile, buf); // nom du domaine
    char* zoneName = new char[BUFSIZE+1];
    strcpy(zoneName, buf);
    zoneNames.push_back(zoneName);
    if (ret == 2) skipLine(ptrFile);
    
    ret = readWord(ptrFile, buf); // type
    if (ret == 2) skipLine(ptrFile);
    if (strcmp(buf, "3D") != 0) // 2D ou autre
      skipLine(ptrFile); // skip definitions 2D
    
    skipLine(ptrFile); // echelle de longueur
    ret = readInt(ptrFile, nvertex);
    if (ret == 1) skipLine(ptrFile);
    //printf("nbre de vertex=%d\n", nvertex);

    ret = readInt(ptrFile, nfaces); 
    if (ret == 1) skipLine(ptrFile);
    //printf("nbre de faces=%d\n", nfaces);

    ret = readInt(ptrFile, ncells); 
    if (ret == 1) skipLine(ptrFile);
    //printf("nbre de cellules=%d\n", ncells);
    
    ret = readInt(ptrFile, nboundaryfaces); 
    if (ret == 1) skipLine(ptrFile);

    //printf("Coords\n");
    FldArrayF* fp = new FldArrayF(nvertex, 3);
    unstructField.push_back(fp);
    E_Float* xp = fp->begin(1); 
    E_Float* yp = fp->begin(2);
    E_Float* zp = fp->begin(3);
    eltType.push_back(8);
    // Coordonnees
    skipLine(ptrFile);
    for (E_Int i = 0; i < nvertex; i++)
    {
      readInt(ptrFile, value);
      readDouble(ptrFile, x); 
      readDouble(ptrFile, y); 
      ret = readDouble(ptrFile, z);
      xp[i] = x; yp[i] = y; zp[i] = z;
      //printf("%f %f %f\n", x,y,z);
    }
    if (ret == 1) skipLine(ptrFile);
    //printf("connect face\n");
    // Connectivite Faces->Noeuds
    // On est oblige de la lire en storage variable
    skipLine(ptrFile);
    E_LONG lpos;
    lpos = KFTELL(ptrFile);
    E_Int size = 0;
    for (E_Int i = 0; i < nfaces; i++)
    {
      readInt(ptrFile, value); // no de la face
      ret = readInt(ptrFile, np); size += np+1;
      for (E_Int j = 0; j < np; j++)
      {
        ret = readInt(ptrFile, value);
      }
    }                 
    KFSEEK(ptrFile, lpos, SEEK_SET);

    //printf("size=%d\n", size);
    FldArrayI connect1(size);
    size = 0;
    
    for (E_Int i = 0; i < nfaces; i++)
    {
      readInt(ptrFile, value); // no de la face
      ret = readInt(ptrFile, np); connect1[size] = np; size++;
      //printf("np=%d\n", np);
      //if (size+np+1 > connect1.getSize()) 
      //  connect1.reAlloc(size+(nfaces-i+10)*(np+1));
      for (E_Int j = 0; j < np; j++)
      {
        ret = readInt(ptrFile, value); connect1[size] = value; size++;
      }
    }
    //connect1.reAlloc(size);
    //printf("%d\n", size);

    // Connectivite faces->elts
    //printf("face elts\n");
    skipLine(ptrFile);
    FldArrayI connect2(nfaces*2);
    E_Int* cn2 = connect2.begin();
    for (E_Int i = 0; i < nfaces; i++)
    {
      readInt(ptrFile, value); //printf("no=%d\n", value); // no face
      readInt(ptrFile, value);  //printf("nv=%d\n", value); // nbre de faces valides (1 ou 2)
      ret = readInt(ptrFile, value); cn2[2*i] = value;
      ret = readInt(ptrFile, value); cn2[2*i+1] = value;
      //printf("%d %d\n", cn2[2*i], cn2[2*i+1]);
    }
    if (ret == 1) skipLine(ptrFile);

    // Construit la connectivite elts->faces
    // Compte les faces pour chaque elements
    //printf("elts face\n");
    FldArrayI nbface(ncells); nbface.setAllValuesAtNull();
    E_Int* nbfacep = nbface.begin();
    E_Int a, b;
    for (E_Int i = 0; i < nfaces; i++)
    {
      a = cn2[2*i]; b = cn2[2*i+1];
      nbfacep[a-1]++;
      if (b != 0) nbfacep[b-1]++;
    }
    //for (E_Int i = 0; i < ncells; i++) printf("%d %d\n", i, nbfacep[i]);

    size = ncells;
    for (E_Int i = 0; i < ncells; i++) size += nbfacep[i];
    //printf("pos %d\n", ncells);
    FldArrayI pos(ncells); // position des faces pour chaque elt
    E_Int* posp = pos.begin();
    posp[0] = 0;
    for (E_Int i = 0; i < ncells-1; i++)
    {
      posp[i+1] = posp[i]+(nbfacep[i]+1);
    }

    //printf("pb %d\n", size);
    FldArrayI connect3(size); E_Int* connect3p = connect3.begin();
    for (E_Int i = 0; i < ncells; i++)
    {
      connect3p[posp[i]] = nbface[i];
    }

    nbface.setAllValuesAtNull();
    for (E_Int i = 0; i < nfaces; i++)
    {
      a = cn2[2*i]; b = cn2[2*i+1];
      connect3p[posp[a-1]+nbfacep[a-1]+1] = i+1; nbfacep[a-1]++;
      if (b != 0) {connect3p[posp[b-1]+nbfacep[b-1]+1] = i+1; nbfacep[b-1]++;}
    }
    pos.malloc(0); connect2.malloc(0); nbface.malloc(0);

    // Fusionne les connectivites
    //printf("fusion\n");
    FldArrayI* cn = new FldArrayI(connect1.getSize()+connect3.getSize()+4, 1);
    connect.push_back(cn);
    E_Int* cp = cn->begin();
    cp[0] = nfaces;
    cp[1] = connect1.getSize();    
    for (E_Int i = 0; i < connect1.getSize(); i++)
    {
      cp[i+2] = connect1[i];
    }
    int pt = 2+connect1.getSize();
    cp[pt] = ncells;
    cp[pt+1] = connect3.getSize();
    for (E_Int i = 0; i < connect3.getSize(); i++)
    {
      cp[pt+2+i] = connect3[i];
    }
    connect3.malloc(0); connect1.malloc(0);

    // Lit les frontieres
    //printf("Lecture frontieres %d\n", nboundaryfaces);
    skipLine(ptrFile);
#define BCSTRINGMAXSIZE 50
    FldArrayI* faces = new FldArrayI(nboundaryfaces);
    char* names = new char [nboundaryfaces * BCSTRINGMAXSIZE]; // taille max des etiquettes de BC
    E_Int* facesp = faces->begin();
    E_Int le;
    E_Int c = 0;

    for (E_Int i = 0; i < nboundaryfaces; i++)
    {
      ret = readInt(ptrFile, value); // numerotation (skip) 
      ret = readInt(ptrFile, value); // no de la face
      facesp[i] = value;
      ret = readWord(ptrFile, buf); // nom de la BC
      le = strlen(buf);
      le = K_FUNC::E_min(le, BCSTRINGMAXSIZE-1);
      for (E_Int l = 0; l < le; l++) { names[c+l] = buf[l]; }
      c += le;
      names[c] = '\0'; c++;
    }
    
    BCFaces.push_back(faces);
    BCNames.push_back(names);
  } // Blocs

  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Write arrays as cedre file.
   Write only NGON_n arrays. */
//=============================================================================
E_Int K_IO::GenIO::cedrewrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames,
  PyObject* BCFaces)
{
  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: cedrewrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

 char format1[40]; char fmtcrd[121]; char dataFmtl[40];
 strcpy(dataFmtl, dataFmt);
 int l = strlen(dataFmt); 
 if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

 // Build format for data
 strcpy(format1,"%d ");
 sprintf(fmtcrd,"%s %s %s\n", dataFmt, dataFmt, dataFmtl);
 strcat(format1,fmtcrd);

 // BCFaces size
 E_Int BCFacesSize = 0;
 if (PyList_Check(BCFaces) == true) BCFacesSize = PyList_Size(BCFaces);
 IMPORTNUMPY;

  // Ecriture de l'entete
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("cedrewrite: I can't open file %s.\n", file);
    return 1;
  }

  // Nbre de domaines
  unsigned int nd = 0;
  for (unsigned int i = 0; i < eltType.size(); i++)
  {
    if (eltType[i] == 8) nd++; 
  }
  if (nd != eltType.size())
    printf("Warning: cedrewrite: array list contain non-NGons arrays. Skipped...\n");

  fprintf(ptrFile, "%d : Nb de dom.\n", nd);

  E_Int eltTypeSize = eltType.size();
  for (E_Int i = 0; i < eltTypeSize; i++)
  {
    if (eltType[i] == 8)
    {
      fprintf(ptrFile, " -----------------\n");
      fprintf(ptrFile, " 1.0 DONNEES GENERALES DU DOMAINE\n");
      fprintf(ptrFile, "   %s\n", zoneNames[i]);
      fprintf(ptrFile, "   3D\n");
      fprintf(ptrFile, "   1.0\n");
      
      FldArrayI& cn = *connect[i];
      FldArrayF& f = *unstructField[i];
      E_Int* cnp = cn.begin();
      E_Int nvertex = f.getSize();
      E_Int nfaces = cnp[0];
      E_Int pt = cnp[1];
      E_Int ncells = cnp[pt+2];
      E_Int facLim = 0;
      if (i < BCFacesSize)
      {
        PyObject* BCs = PyList_GetItem(BCFaces, i);
        E_Int size = PyList_Size(BCs);
        for (E_Int j = 0; j < size/2; j++)
        {
          PyArrayObject* array = (PyArrayObject*)PyList_GetItem(BCs, 2*j+1);
          E_Int np = PyArray_SIZE(array);
          facLim += np;
        }
      }
      fprintf(ptrFile, "     %d       NB NOEUDS\n", nvertex);
      fprintf(ptrFile, "     %d       NB FACES\n", nfaces);
      fprintf(ptrFile, "     %d       NB ELMTS\n", ncells);
      fprintf(ptrFile, "     %d       NB FAC LIM\n", facLim);
      
      // Coordonnees des noeuds
      fprintf(ptrFile, "1. GRID NODES : NODE no., x, y, z\n");
      for (E_Int j = 0; j < nvertex; j++)
      {
        fprintf(ptrFile, format1, 
                j+1, f(j,posx), f(j,posy), f(j,posz)); 
      }

      // Connectivite faces->noeuds
      fprintf(ptrFile, "2. FACES -> NODES : FACE no., number of NODES, no of NODE 1,...\n");
      pt = 2; E_Int nf;
      for (E_Int j = 0; j < nfaces; j++)
      {
        nf = cnp[pt]; pt++;
        fprintf(ptrFile, " %d  %d", j+1, nf);
        for (E_Int k = 0; k < nf; k++)
        { fprintf(ptrFile, " %d", cnp[pt]); pt++; }
        fprintf(ptrFile, "\n");
      }
            
      // Connectivite faces->elts
      FldArrayI cFE;
      K_CONNECT::connectNG2FE(cn, cFE);
      E_Int* facesp1 = cFE.begin(1);
      E_Int* facesp2 = cFE.begin(2);

      fprintf(ptrFile, "3. FACES -> ELTS : FACE no., number of Elts, ELT 1, ELT2\n");
      E_Int jp = 1;
      for (E_Int j = 0; j < nfaces; j++)
      {
        nd = 2;
        if (facesp1[j] == 0 && facesp2[j] == 0)
        {
          nd = 0;
          fprintf(ptrFile, "%d %d %d %d\n", jp, nd, facesp1[j]+1, facesp2[j]+1); jp++; // this is strange!
        }
        else if (facesp1[j] == 0)
        {
          nd = 1;
          fprintf(ptrFile, "%d %d %d %d\n", jp, nd, facesp2[j], facesp1[j]); jp++;
        }
        else if (facesp2[j] == 0)
        {
          nd = 1;
          fprintf(ptrFile, "%d %d %d %d\n", jp, nd, facesp1[j], facesp2[j]); jp++;
        }
        else
        { 
          fprintf(ptrFile, "%d %d %d %d\n", jp, nd, facesp1[j], facesp2[j]); jp++;
        }
      }

      // Faces marquees, a partir de l'objet python BCFaces
      fprintf(ptrFile, "4. FACES MARQUEES : no de face marquee, no de face, name of BC\n");
      
      if (i < BCFacesSize)
      {
        PyObject* BCs = PyList_GetItem(BCFaces, i);
        E_Int size = PyList_Size(BCs);
        E_Int c = 1;
        for (E_Int j = 0; j < size/2; j++)
        {
          char* name = NULL; 
          PyObject* o = PyList_GetItem(BCs, 2*j);
          if (PyString_Check(o)) name = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(o)) name = PyBytes_AsString(PyUnicode_AsUTF8String(o));
#endif
          PyArrayObject* array = (PyArrayObject*)PyList_GetItem(BCs, 2*j+1);
          int* ptr = (int*)PyArray_DATA(array);
          E_Int np = PyArray_SIZE(array);
          for (E_Int k = 0; k < np; k++)
          {
            fprintf(ptrFile, "%d %d %s\n", c, ptr[k], name); c++;
          }
        }
      }
      
    } // NGONS
  }
 
  fclose(ptrFile);
  return 0;
}
