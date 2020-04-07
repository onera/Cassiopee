/*    
    Copyright 2013-2020 Onera.

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

// Binary archive (CEDRE) file support

# include "GenIO.h"
# include <stdio.h>
# include <string.h>
# include "Array/Array.h"
# include <map>

// for dbx
#include <iostream>

using namespace std;
using namespace K_FLD;

class mesh
{

public:
  mesh() { _isomu=NULL; _ifacu=NULL; _icelu=NULL; _pe=NULL; _nbpoints=NULL; 
    _numpoints=NULL; _x=NULL; _y=NULL; _z=NULL; };
  ~mesh() { delete [] _isomu; delete [] _ifacu; delete [] _icelu; 
    delete [] _pe; delete [] _nbpoints; delete [] _numpoints; 
    delete [] _x; delete [] _y; delete [] _z; };
  
public:
  // IN: Numero absolu + temps permette de reperer un bloc de facon unique
  // numero absolu du sous domaine
  unsigned int _nd;
  // temps
  double _temps;
  
  // numero du domaine utilisateur auquel ce bloc appartient
  unsigned int _numutil;
  // nom du du domaine
  char _nom[256];
  // type du maillage
  unsigned char _type;
  // situation des variables (0,1)
  unsigned char _situ;
  // numero du melange
  int _num_mel;
  // numero du groupe
  int _num_grp;
  // nombre de sommets
  int _nsom;
  // nombre de faces
  int _nfac;
  // nombre de cellules limites
  int _nflm;
  int _nflp;
  // nombre de cellules internes
  int _nceli;
  // isomu
  int* _isomu;
  // ifacu
  int* _ifacu;
  // icelu
  int* _icelu;
  // elmin
  unsigned int _elminx;
  // nelem
  unsigned int _nelemx;
  // parent elements
  unsigned int* _pe;
  // nbpoints
  unsigned char* _nbpoints;
  // numpoints
  unsigned int* _numpoints;
  // imax jmax kmax
  unsigned short _imax;
  unsigned short _jmax;
  unsigned short _kmax;
  // elmin
  unsigned short _elmin;
  // nelem
  unsigned short _nelem;
  // comp
  unsigned char _comp_x;
  double _vmin_x;
  double _vmax_x;
  double* _x;
  unsigned char _comp_y;
  double _vmin_y;
  double _vmax_y;
  double* _y;
  unsigned char _comp_z;
  double _vmin_z;
  double _vmax_z;
  double* _z;
  
  };

// Verifie si le domaine nd existe dans meshes, sinon le cree
// il faudrait aussi pouvoir utiliser le temps
mesh* checkMesh(int nd, std::map<unsigned int, mesh*>& meshes)
{
  std::map<unsigned int,mesh*>::iterator it;
  it = meshes.find(nd);
  if (it != meshes.end()) return it->second;
  else
  {
    mesh* m = new mesh; 
    meshes[nd] = m;
    return m;
  }
}

// Read a block name jusqu'au caractere nul
void readBlockName(FILE* ptrFile, char* name)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  printf("blockname: %s\n", name);
}

// Lit la version du fichier
void readVersion(FILE* ptrFile, unsigned char* version)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { version[i] = c; i++; }
  version[i] = '\0';
  printf("version: %u\n", *version);
}

// Lit le titre du fichier
void readTitle(FILE* ptrFile, char* titre, char* date, char* machine, 
               double& tempsZero, double& epsilon)
{
  int c; E_Int i;
  if (titre != NULL)
  {
    i = 0;
    while ((c = fgetc(ptrFile)) != '\0')
    { titre[i] = c; i++; }
    titre[i] = '\0';
  }
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { date[i] = c; i++; }
  date[i] = '\0';
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { machine[i] = c; i++; }
  machine[i] = '\0';
  fread(&tempsZero, sizeof(double), 1, ptrFile); tempsZero = DBE(tempsZero);
  fread(&epsilon, sizeof(double), 1, ptrFile); epsilon = DBE(epsilon);
  
  if (titre != NULL) printf("titre: %s\n", titre);
  printf("date: %s\n", date);
  printf("machine: %s\n", machine);
  printf("temps0: %g\n", tempsZero);
  printf("epsilon: %g\n", epsilon);
}

// Lit les donnes globales
void readGlobal(FILE* ptrFile, unsigned char& solverType, unsigned char& dimField)
{
  fread(&solverType, sizeof(unsigned char), 1, ptrFile);
  fread(&dimField, sizeof(unsigned char), 1, ptrFile);
  printf("solverType: %u\n", solverType);
  printf("dimField: %u\n", dimField);
}

// Lit une espece
void readEspece(FILE* ptrFile, int& numMel, char* nomMel, int& nespeces, char* nomEspece)
{
  int c;
  // numero du melange
  fread(&numMel, sizeof(int), 1, ptrFile); numMel = IBE(numMel);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomMel[i] = c; i++; }
  nomMel[i] = '\0';
  fread(&nespeces, sizeof(int), 1, ptrFile); nespeces = IBE(nespeces);
  i = 0;
  for (E_Int j = 0; j < nespeces; j++)
  {
    while ((c = fgetc(ptrFile)) != '\0')
    { nomEspece[i] = c; i++; }
    nomEspece[i] = '\0';
  }
  printf("numMel: %d\n", numMel);
  printf("nomMel: %s\n", nomMel);
  printf("nbreEspeces: %d\n", nespeces);
  for (E_Int i = 0; i < nespeces; i++) printf("%d %s\n", i, nomEspece);
}

// Lit les donnees scalaires
void readScalar(FILE* ptrFile, int& numGrp, char* nomGrp, int& nelem, 
                char* nomScalar, unsigned char& typeSca)
{
  int c;
  fread(&numGrp, sizeof(int), 1, ptrFile); numGrp = IBE(numGrp);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomGrp[i] = c; i++; }
  nomGrp[i] = '\0';
  fread(&nelem, sizeof(int), 1, ptrFile); nelem = IBE(nelem);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { nomScalar[i] = c; i++; }
  nomScalar[i] = '\0';
  fread(&typeSca, sizeof(unsigned char), 1, ptrFile);
}

// Lit une unite
void readUnit(FILE* ptrFile, char* name, char* unit)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  printf("name %s\n", name);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { unit[i] = c; i++; }
  unit[i] = '\0';
  printf("unit %s\n", unit);
}

// Lit les donnees thermo
void readThermo(FILE* ptrFile, unsigned int& read, char* name,  
  unsigned int& nGe, unsigned int& nGr,
  unsigned int& ne, unsigned int& nr, 
  unsigned int*& dthGe, unsigned int*& dthGr,
  unsigned int*& dthe, double*& dthr)
{
  fread(&read, sizeof(unsigned int), 1, ptrFile); read = UIBE(read);
  printf("thermo: read %u\n", read);
  
  int c;
  E_Int i = 0;
  while ( (c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  printf("thermo name : %s\n", name);
  if (read == 1) return;
  fread(&nGe, sizeof(unsigned int), 1, ptrFile); nGe = UIBE(nGe);
  fread(&nGr, sizeof(unsigned int), 1, ptrFile); nGr = UIBE(nGr);
  fread(&ne, sizeof(unsigned int), 1, ptrFile); ne = UIBE(ne);
  fread(&nr, sizeof(unsigned int), 1, ptrFile); nr = UIBE(nr);
  printf("nGe: %d %d %d %d\n", nGe, nGr, ne, nr);
  dthGe = new unsigned int [nGe];
  dthGr = new unsigned int [nGr];
  dthe = new unsigned int [ne];
  dthr = new double [nr];
  fread(&dthGe, sizeof(unsigned int), nGe, ptrFile); 
  fread(&dthGr, sizeof(unsigned int), nGr, ptrFile);
  fread(&dthe, sizeof(unsigned int), ne, ptrFile);
  fread(&dthr, sizeof(double), nr, ptrFile);
}
 
// Lit la structure d'un maillage
void readStructure(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{  
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  
  fread(&m->_numutil, sizeof(unsigned int), 1, ptrFile); m->_numutil = UIBE(m->_numutil);

  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { m->_nom[i] = c; i++; }
  m->_nom[i] = '\0';
  // type=0 (structure), type=1 (mixte), type=2 (elements finis), type=3 (particules)
  fread(&m->_type, sizeof(unsigned char), 1, ptrFile); 
  fread(&m->_situ, sizeof(unsigned char), 1, ptrFile); 
  fread(&m->_num_mel, sizeof(int), 1, ptrFile); m->_num_mel = IBE(m->_num_mel);
  fread(&m->_num_grp, sizeof(int), 1, ptrFile); m->_num_grp = IBE(m->_num_grp);
  printf("Structure: nd=%u numutil=%u\n", nd, m->_numutil);
  printf("nom dom=%s\n", m->_nom);
  printf("type=%u, situ=%u\n", m->_type, m->_situ);
}

void readDomutil(FILE* ptrFile, unsigned int& numUti, char* nomUti)
{
  int c;
  E_Int i = 0;
  fread(&numUti, sizeof(unsigned int), 1, ptrFile); numUti = UIBE(numUti);
  while ((c = fgetc(ptrFile)) != '\0')
  { nomUti[i] = c; i++; }
  nomUti[i] = '\0';
  printf("nom du domaine %s\n", nomUti);
}

void readNumerotation(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  
  fread(&m->_temps, sizeof(double), 1, ptrFile); m->_temps = DBE(m->_temps);
  
  fread(&m->_nsom, sizeof(int), 1, ptrFile); m->_nsom = IBE(m->_nsom);
  fread(&m->_nfac, sizeof(int), 1, ptrFile); m->_nfac = IBE(m->_nfac);
  fread(&m->_nflm, sizeof(int), 1, ptrFile); m->_nflm = IBE(m->_nflm);
  fread(&m->_nflp, sizeof(int), 1, ptrFile); m->_nflp = IBE(m->_nflp);
  fread(&m->_nceli, sizeof(int), 1, ptrFile); m->_nceli = IBE(m->_nceli);
  printf("%d %d %d\n", m->_nsom, m->_nfac, m->_nflm);
  m->_isomu = new int [m->_nsom];
  fread(m->_isomu, sizeof(int), m->_nsom, ptrFile);
  int* pt = m->_isomu;
  for (E_Int i = 0; i < m->_nsom; i++) pt[i] = IBE(pt[i]);

  m->_ifacu = new int [m->_nfac];
  fread(m->_ifacu, sizeof(int), m->_nfac, ptrFile);
  pt = m->_ifacu;
  for (E_Int i = 0; i < m->_nfac; i++) pt[i] = IBE(pt[i]);
  
  E_Int size = m->_nceli+m->_nflm+m->_nflp;
  m->_icelu = new int [size];
  fread(m->_icelu, sizeof(int), size, ptrFile);
  pt = m->_icelu;
  for (E_Int i = 0; i < size; i++) pt[i] = IBE(pt[i]);
}

void readVariables(FILE* ptrFile, unsigned int& numUti, unsigned int& nvar, char* nomVar,
  char* nomGrp, char* nomCat, unsigned char& glob, unsigned int& elmin,
  unsigned int& nelem, double*& data)
{
  fread(&numUti, sizeof(unsigned int), 1, ptrFile); numUti = UIBE(numUti);
  fread(&nvar, sizeof(unsigned int), 1, ptrFile); nvar = UIBE(nvar);
  
}

void readConnexion(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  
  fread(&m->_temps, sizeof(double), 1, ptrFile); m->_temps = DBE(m->_temps);
  fread(&m->_elminx, sizeof(unsigned int), 1, ptrFile); m->_elminx = UIBE(m->_elminx);
  fread(&m->_nelemx, sizeof(unsigned int), 1, ptrFile); m->_nelemx = UIBE(m->_nelemx);
  
  E_Int size = 2*m->_nelemx;
  m->_pe = new unsigned int [size];
  fread(m->_pe, sizeof(unsigned int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) m->_pe[i] = UIBE(m->_pe[i]);
  
  m->_nbpoints = new unsigned char [m->_nelemx];
  fread(m->_nbpoints, sizeof(unsigned char), m->_nelemx, ptrFile);
  
  E_Int nbPointsTot = 0;
  for (size_t i = 0; i < m->_nelemx; i++) nbPointsTot += (E_Int)m->_nbpoints[i];
    
  size = nbPointsTot;
  m->_numpoints = new unsigned int [size];
  fread(m->_numpoints, sizeof(unsigned int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) m->_numpoints[i] = UIBE(m->_numpoints[i]);
}

// Lit un maillage structure (type=0 ou type=1)
void readMaillage0(FILE* ptrFile, mesh* m)
{
  fread(&m->_temps, sizeof(double), 1, ptrFile); m->_temps = DBE(m->_temps);
  fread(&m->_imax, sizeof(unsigned short), 1, ptrFile); 
  fread(&m->_jmax, sizeof(unsigned short), 1, ptrFile); 
  fread(&m->_kmax, sizeof(unsigned short), 1, ptrFile);
  fread(&m->_comp_x, sizeof(unsigned short), 1, ptrFile); 
  if (m->_comp_x != 0)
  {
    fread(&m->_vmin_x, sizeof(double), 1, ptrFile); m->_vmin_x = DBE(m->_vmin_x);
  }
  if (m->_comp_x != 0 && m->_comp_x != 4)
  {
    fread(&m->_vmax_x, sizeof(double), 1, ptrFile); m->_vmax_x = DBE(m->_vmax_x);
  }
  m->_x = new double [m->_imax*m->_jmax*m->_kmax];
  fread(m->_x, sizeof(double), m->_imax*m->_jmax*m->_kmax, ptrFile);
  for (E_Int i=0; i < m->_imax*m->_jmax*m->_kmax; i++) m->_x[i] = DBE(m->_x[i]); 

  fread(&m->_comp_y, sizeof(unsigned short), 1, ptrFile); 
  if (m->_comp_y != 0)
  {
    fread(&m->_vmin_y, sizeof(double), 1, ptrFile); m->_vmin_y = DBE(m->_vmin_y);
  }
  if (m->_comp_y != 0 && m->_comp_y != 4)
  {
    fread(&m->_vmax_y, sizeof(double), 1, ptrFile); m->_vmax_y = DBE(m->_vmax_y);
  }
  m->_y = new double [m->_imax*m->_jmax*m->_kmax];
  fread(m->_y, sizeof(double), m->_imax*m->_jmax*m->_kmax, ptrFile);
  for (E_Int i=0; i < m->_imax*m->_jmax*m->_kmax; i++) m->_y[i] = DBE(m->_y[i]); 

  fread(&m->_comp_z, sizeof(unsigned short), 1, ptrFile);
  if (m->_comp_z != 0)
  {
    fread(&m->_vmin_z, sizeof(double), 1, ptrFile); m->_vmin_z = DBE(m->_vmin_z);
  }
  if (m->_comp_z != 0 && m->_comp_z != 4)
  {
    fread(&m->_vmax_z, sizeof(double), 1, ptrFile); m->_vmax_z = DBE(m->_vmax_z);
  }
  m->_z = new double [m->_imax*m->_jmax*m->_kmax];
  fread(m->_z, sizeof(double), m->_imax*m->_jmax*m->_kmax, ptrFile);
  for (E_Int i=0; i < m->_imax*m->_jmax*m->_kmax; i++) m->_z[i] = DBE(m->_z[i]); 
}

// Lit un maillage non structure (type=2 ou type=3)
void readMaillage1(FILE* ptrFile, mesh* m)
{
  fread(&m->_temps, sizeof(double), 1, ptrFile); m->_temps = DBE(m->_temps);
  printf("temps=%f\n", m->_temps);
  fread(&m->_elmin, sizeof(unsigned short), 1, ptrFile); 
  fread(&m->_nelem, sizeof(unsigned short), 1, ptrFile);
  fread(&m->_comp_x, sizeof(unsigned short), 1, ptrFile);
  if (m->_comp_x != 0)
  {
    fread(&m->_vmin_x, sizeof(double), 1, ptrFile); m->_vmin_x = DBE(m->_vmin_x);
  }
  if (m->_comp_x != 0 && m->_comp_x != 4)
  {
    fread(&m->_vmax_x, sizeof(double), 1, ptrFile); m->_vmax_x = DBE(m->_vmax_x);
  }
  m->_x = new double [m->_nelem];
  fread(m->_x, sizeof(double), m->_nelem, ptrFile);
  for (E_Int i=0; i < m->_nelem; i++) m->_x[i] = DBE(m->_x[i]); 
    
  fread(&m->_comp_y, sizeof(unsigned short), 1, ptrFile); 
  if (m->_comp_y != 0)
  {
    fread(&m->_vmin_y, sizeof(double), 1, ptrFile); m->_vmin_y = DBE(m->_vmin_y);
  }
  if (m->_comp_y != 0 && m->_comp_y != 4)
  {
    fread(&m->_vmax_y, sizeof(double), 1, ptrFile); m->_vmax_y = DBE(m->_vmax_y);
  }
  m->_y = new double [m->_nelem];
  fread(m->_y, sizeof(double), m->_nelem, ptrFile);
  for (E_Int i=0; i < m->_nelem; i++) m->_y[i] = DBE(m->_y[i]); 

  fread(&m->_comp_z, sizeof(unsigned short), 1, ptrFile); 
  if (m->_comp_z != 0)
  {
    fread(&m->_vmin_z, sizeof(double), 1, ptrFile); m->_vmin_z = DBE(m->_vmin_z);
  }
  if (m->_comp_z != 0 && m->_comp_z != 4)
  {
    fread(&m->_vmax_z, sizeof(double), 1, ptrFile); m->_vmax_z = DBE(m->_vmax_z);
  }
  m->_z = new double [m->_nelem];
  fread(m->_z, sizeof(double), m->_nelem, ptrFile);
  for (E_Int i=0; i < m->_nelem; i++) m->_z[i] = DBE(m->_z[i]); 
  
}

// Lecture generale de maillage
void readMaillage(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  printf("reading type=%d for mesh no %d\n", m->_type, nd);
  if (m->_type == 0 || m->_type == 1) readMaillage0(ptrFile, m);
  else readMaillage1(ptrFile, m);
}

//=============================================================================
/* 
   arcread - read binary cedre format (archive)
   IN: file: file name,
   OUT: varString: variables string
   OUT: structField: field for each structured zones,
   OUT: ni, nj, nk: number of points of each structured zones,
   OUT: unstructField: field for each unstructured zones, 
   OUT: connectivity: connectivity for each unstructured zones,
   OUT: eltType: eltType for each unstructured zones.

   eltType is:
   1: BAR
   2: TRI
   3: QUAD
   4: TETRA
   7: HEXA
   8: NGON
   return 1 if failure.
   return 0 if ok.

*/
//=============================================================================
E_Int K_IO::GenIO::arcread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connectivity,
  vector<E_Int>& eltType, vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: arcread: cannot open file %s.\n", file);
    return 1;
  }

  //printf("Error: arcread: not implemented.\n");
  unsigned char version;
  char name[256];
  char fileName[256];
  char unitName[256];
  char nomUti[256];
  char titre[256];
  char date[256];
  char machine[256];
  double tempsZero;
  double epsilon;
  unsigned char solverType;
  unsigned char dimField;
  int numMel;
  char nomMel[256];
  int nespeces;
  char nomEspece[256*10];
  int numGrp;
  char nomGrp[256];
  int nelem;
  char nomScalar[256];
  unsigned char typeSca;
  char unit[256];
  unsigned int read;
  unsigned int numUti;
  unsigned int nGe, nGr, ne, nr;
  unsigned int* dthGe=NULL;
  unsigned int *dthGr=NULL;
  unsigned int *dthe=NULL;
  double *dthr=NULL;
  int* isomu=NULL;
  int* ifacu=NULL;
  int* icelu=NULL;
  unsigned int* num=NULL;
  unsigned char* nbpoints=NULL;
  unsigned int* numpoints=NULL;
  std::map<unsigned int, mesh*> meshes;
  
  readBlockName(ptrFile, name);
  readVersion(ptrFile, &version);
  
  while (! feof(ptrFile))
  {
    readBlockName(ptrFile, name);
    if (strcmp(name, "VERSION") == 0) readVersion(ptrFile, &version);
    else if (strcmp(name, "TITRE") == 0) readTitle(ptrFile, titre, date, machine, tempsZero, epsilon);
    else if (strcmp(name, "SANS TITRE") == 0) readTitle(ptrFile, NULL, date, machine, tempsZero, epsilon);
    else if (strcmp(name, "GLOBAL") == 0) readGlobal(ptrFile, solverType, dimField);
    else if (strcmp(name, "ESPECES") == 0) readEspece(ptrFile, numMel, nomMel, nespeces, nomEspece);
    else if (strcmp(name, "SCALAIRES") == 0) readScalar(ptrFile, numGrp, nomGrp, nelem, nomScalar, typeSca);
    else if (strcmp(name, "UNITE") == 0) readUnit(ptrFile, unitName, unit);
    else if (strcmp(name, "THERMODYNAMIQUE") == 0) readThermo(ptrFile, read, fileName, nGe, nGr, ne, nr, dthGe, dthGr, dthe, dthr);
    else if (strcmp(name, "DOMUTIL") == 0) readDomutil(ptrFile, numUti, nomUti);
    else if (strcmp(name, "STRUCTURE") == 0) readStructure(ptrFile, meshes);
    else if (strcmp(name, "NUMEROTATION") == 0) readNumerotation(ptrFile, meshes);
    else if (strcmp(name, "CONNEXION") == 0) readConnexion(ptrFile, meshes);
    else if (strcmp(name, "MAILLAGE") == 0) readMaillage(ptrFile, meshes);
    else { printf("Block pas encore implemente: %s\n", name); break; }
  }
  delete [] dthGe; delete [] dthGr; delete [] dthe; delete [] dthr;
  delete [] isomu; delete [] ifacu; delete [] icelu;
  delete [] num; delete [] nbpoints; delete [] numpoints;
  return 0;
}

//=============================================================================
/*
  This routine enables binary write of archive file.
*/
//=============================================================================
E_Int K_IO::GenIO::arcwrite(char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector<E_Int>& eltType,
      std::vector<char*>& zoneNames)
{
  printf("Error: arcwrite: not implemented.\n");
  return 1;
}

