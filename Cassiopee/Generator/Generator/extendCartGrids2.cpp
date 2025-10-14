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

#include "generator.h"

using namespace K_FLD;
using namespace std;
using namespace K_CONST;

//=============================================================================
/* Extend Cartesian grids w.r.t depth with a minimum overlapping */
//=============================================================================
PyObject* K_GENERATOR::extendCartGrids2(PyObject* self, PyObject* args)
{
  PyObject *arrays;
  E_Int ext, optimized, extBnd;
  if (!PYPARSETUPLE_(args, O_ III_, &arrays, &ext, &optimized, &extBnd)) return NULL;

  if (ext < 0) 
  {
   PyErr_SetString(PyExc_TypeError, 
                   "extendCartGrids: ext must be a positive value.");
   return NULL;
  }
  if (ext == 0) return arrays;
  if (optimized != 0 && optimized != 1 && optimized != -1) 
  { printf("Warning: extendCartGrids: optimized is set to 1.\n"); optimized = 1; }

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  if ( isOk == -1 ) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extendCartGrids: invalid list of arrays.");
     return NULL;   
  }
  /* Verification des positions de x,y,z */
  E_Int nzones = structF.size();
  E_Int posxi, posyi, poszi;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  for (E_Int i = 0; i < nzones; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]);
    if ( posxi == -1 || posyi == -1 || poszi == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "extendCartGrids: arrays must contain coordinates.");
      for (E_Int v = 0 ; v < nzones; v++) RELEASESHAREDS(objst[v], structF[v]);
      return NULL;
    }
    posxi++; posyi++; poszi++;
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi);
  }
  E_Int dim = 3;
  if ( nkt[0] == 1 ) dim = 2;

  //printf("extend: %d \n", dim);

  // determination des elts dont les bbox intersectent l elt courant
  FldArrayF bbox(nzones, 6);// xmin, ymin, zmin, xmax, ymax, zmax
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  E_Float tol = 1.e-6; E_Float tol2 = tol*tol;
  for (E_Int v = 0; v < nzones; v++)
  {
    K_COMPGEOM::boundingBoxStruct(nit[v], njt[v], nkt[v],
                                  structF[v]->begin(posxt[v]),
                                  structF[v]->begin(posyt[v]), 
                                  structF[v]->begin(poszt[v]),
                                  xminp[v], yminp[v], zminp[v],
                                  xmaxp[v], ymaxp[v], zmaxp[v]);

    minB[0] = xminp[v]; minB[1] = yminp[v]; minB[2] = zminp[v];
    maxB[0] = xmaxp[v]; maxB[1] = ymaxp[v]; maxB[2] = zmaxp[v];
    if (dim == 2) { minB[2] = 0.; maxB[2] = 1.; zminp[v]=0; zmaxp[v]=1; }
    //printf("ymin %f  ymax %f , Noz: %d \n",  yminp[v],ymaxp[v], v );
  }

  // Determination des extensions pour chq zone a partir de l'octree
  E_Int extg = ext; E_Int extf = ext; E_Int extff = ext;
  if (optimized == 1) {extg = ext-1; extff = extg;}
  if (optimized ==-1) {extg = ext+1;}
  //if (optimized ==-1) {extg = ext+1;extff = extg;}

  E_Int mx_sf=20; //nombre maximal de sous fnetre pour raccord d'une face
  FldArrayI voisin(nzones*6*mx_sf);
  FldArrayF deltax(nzones); E_Float* ipt_dx = deltax.begin(1);
  FldArrayI  sousface(nzones*6); sousface.setAllValuesAtNull();
  FldArrayI  resolution(nzones*6); resolution.setAllValuesAtNull();
  FldArrayI extension(nzones, 6); extension.setAllValuesAtNull();
  vector<E_Int> indicesBB;
  E_Int indA1, indB1, indC1, indD1, indE1, indF1, indG1, indH1;
  E_Int indA2, indB2, indC2, indD2, indE2, indF2, indG2, indH2, v2;
  E_Int* ext1 = extension.begin(1);
  E_Int* ext2 = extension.begin(2);
  E_Int* ext3 = extension.begin(3);
  E_Int* ext4 = extension.begin(4);
  E_Int* ext5 = extension.begin(5);
  E_Int* ext6 = extension.begin(6);

  E_Int* iptvois= voisin.begin(1);
  E_Int*  sf    = sousface.begin(1);
  E_Int* resol  = resolution.begin(1);

  E_Int ret;
  E_Float xp, yp, zp, dx, dy, dz;
  E_Float s1, s2; // surface de la facette initiale et surface de facette projetee
  E_Int nbboxes;
  E_Int ni1, nj1, nk1, ni2, nj2, nk2;
  E_Int shift1, shift2;
  E_Float *xt1, *yt1, *zt1, *xt2, *yt2, *zt2;
  E_Float dhmax, dhmin;
  E_Float dymin, dymax, dxmin, dxmax;
  E_Float  xmax2, xmin2, xmax1, xmin1;
  E_Float  ymax2, ymin2, ymax1, ymin1;
  E_Float  zmax2, zmin2, zmax1, zmin1;
  E_Int found1, found2, found3, found4;
  vector< vector<E_Int> > dejaVu(nzones);
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
  E_Float diff, diff1;
  if (dim == 2) 
  {
    FldArrayF face(2,3);//3D
    FldArrayI cnf(1,2); // cn BAR de la facette
    cnf(0,1) = 1; cnf(0,2) = 2;  


    //FldArrayF face(4,3);//facette = 2 TRI
    //FldArrayI cnf(2,3);// cn TRI de la facette de l HEXA
    cnf(0,1) = 1; cnf(0,2) = 2;// cnf(0,3) = 4; cnf(1,1) = 2; cnf(1,2) = 3; cnf(1,3) = 4;

    //E_Int nzones_loc=1;
    E_Int nzones_loc=nzones;
    for (E_Int v1 = 0; v1 < nzones_loc; v1++)
    {
    printf(" properties 2D %d \n", v1);
#include "faceproperties2D.cpp"
    }

    //on etend d'abord les faces possédant plusieurs sous faces (sf1)
    for (E_Int v1 = 0; v1 < nzones_loc; v1++)
    {
     for (E_Int idir = 0; idir < 4; idir++)
     {
       E_Int idir_opp= idir+1;
       if (idir==1 || idir==3 || idir==5) idir_opp= idir-1;
       if (sf[v1*6+ idir] >1)
        {
         
    printf(" multi 2D: v1= %d  ext= %d \n", v1, ext);
           E_Int fullMono=99;
           if (resol[v1*6+idir] ==2 ||  resol[v1*6+idir]==3) // souface grossiere  ou fin et grossiere
             {
               fullMono=1;
               E_Int ext_loc=ext;
               for (E_Int ssf = 0; ssf < sf[v1*6+ idir]; ssf++) //loop sousface
                {
                   E_Int v2 = iptvois[v1*6*mx_sf+ idir*mx_sf + ssf];  //no zone de raccord de la sousface
                   if (sf[v2*6+ idir_opp] >1) fullMono=0;
                   // voisin grossier + 
                   //if (ipt_dx[v2]/ipt_dx[v1] >1.1 && ext1[v2 + idir_opp*nzones] > ext_loc && sf[v2+ idir_opp] >1 ) ext_loc = ext1[v2 + idir_opp*nzones];
                   if (ipt_dx[v2]/ipt_dx[v1] >1.1  ) ext_loc = ext+1;
                }

               ext_loc = ext;
               if (fullMono==1) ext1[v1 +idir*nzones]  = ext_loc;
               else
                 { ext1[v1+idir*nzones]  = ext_loc;
                 }

             }
           else{
                ext1[v1 + idir*nzones]  = ext;
               }
           printf("Nd: %d , idir: %d, sous face: %d , Resol: %d ghost: %d , fullmono: %d \n", v1, idir, sf[v1*6+ idir], resol[v1*6+idir], ext1[v1 + idir*nzones], fullMono);
        }//test sous face
     } //loop idir
    } //loop zone


    //on etend ensuite les faces possédant une seule sous faces (sf1)
    for (E_Int v1 = 0; v1 < nzones_loc; v1++)
    {
     for (E_Int idir = 0; idir < 4; idir++)
     {
       E_Int idir_opp= idir+1;
       if (idir==1 || idir==3 || idir==5) idir_opp= idir-1;
       if (sf[v1*6+ idir] ==1)
        {
            printf(" mono %d \n", v1);
            E_Int v2 = iptvois[v1*6*mx_sf+ idir*mx_sf + 0];
            E_Int ext_loc=ext;
            //printf("v12 %d %d idir: %d , idir_opp: %d , adr: %d , nzones: %d \n", v1, v2, idir, idir_opp,  v1*6*mx_sf+ idir_opp*mx_sf,  nzones  );
             
            if (ipt_dx[v2]/ipt_dx[v1] > 1.1 )
                { 
                 if      (ext1[v2 + idir_opp*nzones]<=ext)   ext_loc +=1;
                 else if (ext1[v2 + idir_opp*nzones]==ext+1) ext_loc -=1;
                  
                  printf("ext_loc %d , idir: %d ext_opp: %d \n", ext_loc, idir, ext1[v2 + idir_opp*nzones]);
                }
            ext_loc = ext;
            ext1[v1+idir*nzones]  = ext_loc;
            printf("Nd: %d , idir: %d, sous face: %d , Resol: %d ghost: %d \n", v1, idir, sf[v1*6+ idir], resol[v1*6+idir], ext1[v1+ idir*nzones]);
            //printf("Nd: %d , idir: %d, sous face: %d , Resol: %d ghost: %d \n", v1, idir, sf[v1*6+ idir], resol[v1*6+idir], ext);
        }
    }
   }//loop zone


  }
  else //( dim == 3 )
  {    
    FldArrayF face(4,3);//facette = 2 TRI
    FldArrayI cnf(2,3);// cn TRI de la facette de l HEXA
    cnf(0,1) = 1; cnf(0,2) = 2; cnf(0,3) = 4;
    cnf(1,1) = 2; cnf(1,2) = 3; cnf(1,3) = 4;

    //E_Int nzones_loc=1;
    E_Int nzones_loc=nzones;
    for (E_Int v1 = 0; v1 < nzones_loc; v1++)
    {
    printf(" properties %d \n", v1);
#include "faceproperties3D.cpp"
    }
    //on etend d'abord les faces possédant plusieurs sous faces (sf1)
    for (E_Int v1 = 0; v1 < nzones_loc; v1++)
    {
     for (E_Int idir = 0; idir < 6; idir++)
     {
       E_Int idir_opp= idir+1;
       if (idir==1 || idir==3 || idir==5) idir_opp= idir-1;
       if (sf[v1*6+ idir] >1)
        {
         
    printf(" multi %d \n", v1);
           E_Int fullMono=99;
           if (resol[v1*6+idir] ==2 ||  resol[v1*6+idir]==3) // souface grossiere  ou fin et grossiere
             {
               fullMono=1;
               E_Int ext_loc=ext;
               for (E_Int ssf = 0; ssf < sf[v1*6+ idir]; ssf++) //loop sousface
                {
                   E_Int v2 = iptvois[v1*6*mx_sf+ idir*mx_sf + ssf];  //no zone de raccord de la sousface
                   if (sf[v2*6+ idir_opp] >1) fullMono=0;
                   // voisin grossier + 
                   //if (ipt_dx[v2]/ipt_dx[v1] >1.1 && ext1[v2 + idir_opp*nzones] > ext_loc && sf[v2+ idir_opp] >1 ) ext_loc = ext1[v2 + idir_opp*nzones];
                   if (ipt_dx[v2]/ipt_dx[v1] >1.1  ) ext_loc = ext+1;
                }

            ext_loc = ext;
               //if (fullMono==1) ext1[v1 +idir*nzones]  = ext;
               if (fullMono==1) ext1[v1 +idir*nzones]  = ext_loc;
               else
                 { ext1[v1+idir*nzones]  = ext_loc;
                 }

             }
           //else if (resol[v1*6+idir] ==-2 ) // souface fine
           //  {
           //    if (fullMono==1) ext1[v1 +idir*nzones]  = ext;
           //  }
           else{
                ext1[v1 + idir*nzones]  = ext;
               }
           printf("Nd: %d , idir: %d, sous face: %d , Resol: %d ghost: %d , fullmono: %d \n", v1, idir, sf[v1*6+ idir], resol[v1*6+idir], ext1[v1 + idir*nzones], fullMono);
        }//test sous face
     } //loop idir
    } //loop zone

    //on etend ensuite les faces possédant une seule sous faces (sf1)
    for (E_Int v1 = 0; v1 < nzones_loc; v1++)
    {
     for (E_Int idir = 0; idir < 6; idir++)
     {
       E_Int idir_opp= idir+1;
       if (idir==1 || idir==3 || idir==5) idir_opp= idir-1;
       if (sf[v1*6+ idir] ==1)
        {
            printf(" mono %d \n", v1);
            E_Int v2 = iptvois[v1*6*mx_sf+ idir*mx_sf + 0];
            E_Int ext_loc=ext;
            //printf("v12 %d %d idir: %d , idir_opp: %d , adr: %d , nzones: %d \n", v1, v2, idir, idir_opp,  v1*6*mx_sf+ idir_opp*mx_sf,  nzones  );
             
            if (ipt_dx[v2]/ipt_dx[v1] > 1.1 )
                { 
                 if      (ext1[v2 + idir_opp*nzones]<=ext)   ext_loc +=1;
                 else if (ext1[v2 + idir_opp*nzones]==ext+1) ext_loc -=1;
                  
                  printf("ext_loc %d , idir: %d ext_opp: %d \n", ext_loc, idir, ext1[v2 + idir_opp*nzones]);
                }
          /*
            else if (ipt_dx[v2]/ipt_dx[v1] < 1.1  && ipt_dx[v2]/ipt_dx[v1] > 0.9)// iso resolution
               { 
                 if (ext1[v2 + idir_opp*nzones]==ext+1) ext_loc -=1;
               }
            */
            ext_loc = ext;
            ext1[v1+idir*nzones]  = ext_loc;
            printf("Nd: %d , idir: %d, sous face: %d , Resol: %d ghost: %d \n", v1, idir, sf[v1*6+ idir], resol[v1*6+idir], ext1[v1+ idir*nzones]);
            //printf("Nd: %d , idir: %d, sous face: %d , Resol: %d ghost: %d \n", v1, idir, sf[v1*6+ idir], resol[v1*6+idir], ext);
        }
    }
   }//loop zone

  }//2d/3d

  PyObject* l = PyList_New(0); 

  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int nk = nkt[v]; 
    E_Float* xp = structF[v]->begin(posxt[v]);
    E_Float* yp = structF[v]->begin(posyt[v]);
    E_Float* zp = structF[v]->begin(poszt[v]);
    E_Int nfldo = structF[v]->getNfld();
    E_Float eps_local = 1.0e-12;
    E_Float dh  = xp[1]-xp[0];
    E_Float dh2 = yp[ni]-yp[0];
    E_Float dh3 = dh;
    if (dim == 3) dh3 = zp[ni*nj]-zp[0];

    //Needed to guarantee the same indices in the tc (pointlist, pointlistdonor, etc.)
    //when dh2 and dh3 are almost the same as dh. E.g. pointlist will be different when
    //dh-dh2=~ 1e-16
    if (abs(dh-dh2)<eps_local) dh2=dh;
    if (abs(dh-dh3)<eps_local) dh3=dh;
    

    if (extBnd > 0) 
    {
      if ( ext1[v] == 0 && extBnd>0) ext1[v]=extBnd;
      if ( ext2[v] == 0 && extBnd>0) ext2[v]=extBnd;
      if ( ext3[v] == 0 && extBnd>0) ext3[v]=extBnd;
      if ( ext4[v] == 0 && extBnd>0) ext4[v]=extBnd;
      if ( ext5[v] == 0 && extBnd>0) ext5[v]=extBnd;
      if ( ext6[v] == 0 && extBnd>0) ext6[v]=extBnd;
    }

    E_Float xxor = xp[0]-ext1[v]*dh;
    E_Float yyor = yp[0]-ext3[v]*dh2;
    E_Float zzor = zp[0]-ext5[v]*dh3;
    RELEASESHAREDS(objst[v], structF[v]);
    E_Int nio = ni+ext1[v]+ext2[v]; E_Int njo = nj+ext3[v]+ext4[v]; E_Int nko = nk+ext5[v]+ext6[v];
    E_Int npts = nio*njo*nko;
    E_Int api = 1;//api 2 plante
    PyObject* tpl = K_ARRAY::buildArray2(nfldo, structVarString[v], nio, njo, nko, api); 
    E_Float* fptr = K_ARRAY::getFieldPtr(tpl);
    FldArrayF newcoords(npts,nfldo, fptr, true);
    E_Float* xn = newcoords.begin(1);
    E_Float* yn = newcoords.begin(2);
    E_Float* zn = newcoords.begin(3);
    E_Int nionjo = nio*njo;
    for (E_Int k = 0; k < nko; k++)    
      for (E_Int j = 0; j < njo; j++)
        for (E_Int i = 0; i < nio; i++)
        {
          E_Int ind = i + j*nio + k*nionjo; 
          xn[ind] = xxor + i*dh;
          yn[ind] = yyor + j*dh2;
          zn[ind] = zzor + k*dh3;
        }
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
      PyObject* extentN = K_NUMPY::buildNumpyArray(extension,1);
  PyObject* tupleOut = Py_BuildValue("[OO]", l, extentN);
      //return l;
  Py_DECREF(l); 
  Py_DECREF(extentN); 
  return tupleOut;
}

