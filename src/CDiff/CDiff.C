//=============================================================================
// This is the differ of output.plt files
// (c) -ONERA- in 2004-2007
//=============================================================================
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream.h>

#include "Gen/IO/GenIO.h"
#include <vector.h>
#include "Fld/Base/FldArray.h"
#include "Blk/Base/BlkMesh.h"
#include "Blk/Interp/BlkInterpAdt.h"

#ifdef __MACH__
// For Mac OSX
#define isnan __isnand
#define isinf __isinfd
#endif

using namespace E_FUNC;

void diff1(E_Int argc, char** argv);
void diff2(E_Int argc, char** argv);

//=============================================================================
// diff1
//=============================================================================
void diff1(E_Int argc, char** argv)
{
  char varString[256];
  char filemr[256];
  char filesr[256];
  char filesm[256];
  strcpy(filemr, argv[1]);
  strcpy(filesr, argv[2]);
  strcpy(filesm, argv[3]);
  
  /*------------------------------------------*/
  /* Reading filemr : reference mesh in nodes */
  /*------------------------------------------*/
  vector<FldArrayF*> fieldmr;
  vector<E_Int> immr;
  vector<E_Int> jmmr;
  vector<E_Int> kmmr;
  vector<FldArrayI*> cmr;
  vector<FldArrayF*> ufieldmr;
  vector<E_Int> etmr;
  cout << "INFO : << Reading reference mesh in nodes : "
       <<filemr<<"..."<<flush;
  GenIO::getInstance()->tecread(filemr, varString, fieldmr, immr, jmmr, kmmr,
                                ufieldmr, cmr, etmr);
  cout << "done."<<endl;
  
  /*------------------------------------------------*/
  /* Reading filesr : reference solution in centers */
  /*------------------------------------------------*/
  vector<FldArrayF*> fieldsr;
  vector<E_Int> imsr;
  vector<E_Int> jmsr;
  vector<E_Int> kmsr;
  vector<FldArrayI*> csr;
  vector<FldArrayF*> ufieldsr;
  vector<E_Int> etsr;
  cout << "INFO : << Reading reference solution in centers : "
       <<filesr<<"..."<<flush;
  GenIO::getInstance()->tecread(filesr, varString, fieldsr, imsr, jmsr, kmsr, 
                                ufieldsr, csr, etsr);
  cout << "done."<<endl;
  
  /*-----------------------------------------------------*/
  /* Reading filesm : reading tested solution in centers */
  /*-----------------------------------------------------*/
  vector<FldArrayF*> fieldsm;
  vector<E_Int> imsm;
  vector<E_Int> jmsm;
  vector<E_Int> kmsm;
  vector<FldArrayI*> csm;
  vector<FldArrayF*> ufieldsm;
  vector<E_Int> etsm;

  cout << "INFO : << Reading tested solution in centers : "
       <<filesm<<"..."<<flush;
  GenIO::getInstance()->tecread(filesm, varString, fieldsm, imsm, jmsm, kmsm, 
                                ufieldsm, csm, etsm);
  cout << "done."<<endl;
 
  /*------------------------------------*/
  /* Building the precond for all grids */
  /*------------------------------------*/
  cout << "INFO : Preconditionning..."<<flush;
  vector<BlkInterpAdt*> adts;
  vector<FldArrayF*> errors;
  E_Int nzone = 0;
  for (vector<FldArrayF*>::iterator it = fieldmr.begin();
       it != fieldmr.end(); it++)
  {
    FldArrayF* f = *it;
    E_Int ni = immr[nzone];
    E_Int nj = jmmr[nzone];
    E_Int nk = kmmr[nzone];
    FldNodeF fn(ni*nj*nk,3);
    for (E_Int ind = 0; ind < ni*nj*nk; ind++)
    {
      fn(ind, 1) = (*f)(ind, 1);
      fn(ind, 2) = (*f)(ind, 2);
      fn(ind, 3) = (*f)(ind, 3);
    }
      
    BlkMesh mesh(ni, nj, nk, fn);
    BlkMesh* meshc = new(BlkMesh);
    meshc->createExtendedCenterMesh(mesh);
    BlkInterpAdt* myAdt = new BlkInterpAdt(*meshc);
    adts.push_back(myAdt);
    nzone++;
  }
  cout << "done."<<endl;
  
  /*----------------------------------------*/
  /* Interpolating from one grid to another */
  /*----------------------------------------*/
  /* Beware : performed interpolation are order 2 so add an error
     to your solution */
  cout << "INFO : Interpolating..."<<flush;
  E_Int ic, jc, kc;
  FldArrayF cf(8);
  E_Int alpha, beta, gamma;
  E_Int ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7;
    
  for (vector<FldArrayF*>::iterator it = fieldsm.begin();
       it != fieldsm.end(); it++)
  {
    FldArrayF& f = **it;
    FldArrayF* error = new FldArrayF(f.getSize(), f.getNfld());
      
    error->setAllValuesAt(1.e6);
    
    for (E_Int ind = 0; ind < f.getSize(); ind++)
    {
      E_Float x = f(ind,1);
      E_Float y = f(ind,2);
      E_Float z = f(ind,3);
      
      (*error)(ind,1) = x;
      (*error)(ind,2) = y;
      (*error)(ind,3) = z;

      E_Int no = 0;
      for (vector<BlkInterpAdt*>::iterator it2 = adts.begin();
           it2 != adts.end(); it2++)
      {
        BlkInterpAdt* adt = *it2;
        E_Boolean found = 
          adt->searchInterpolationCell(x, y, z, ic, jc, kc, cf);
        FldArrayF& f2 = *fieldsr[no];
        E_Int imr = imsr[no];
        E_Int jmr = jmsr[no];
        E_Int kmr = kmsr[no];
        E_Int nic = imr+2;
        E_Int njc = jmr+2;
        E_Int nkc = kmr+2;
        
        if (found == E_True)
        {
          alpha = 1;
          beta = 1;
          gamma = 1;
            
          ic = ic-1;
          jc = jc-1;
          kc = kc-1;
            
          if (ic == 0)
          {
            ic = 1;
            alpha = 0;
          }
          else if (ic == nic-2)
            alpha = 0;
          
          if (jc == 0)
          {
            jc = 1;
            beta = 0;
          }
          else if (jc == njc-2)
            beta = 0;
          
          if (kc == 0)
          {
            kc = 1;
            gamma = 0;
          }
          else if (kc == nkc-2)
            gamma = 0;
            
          // Index of interpolation cell (in center)
          ind0 = ic-1+(jc-1)*imr+(kc-1)*imr*jmr;
          ind1 = ic+alpha-1+(jc-1)*imr+(kc-1)*imr*jmr;
          ind2 = ic-1+(jc+beta-1)*imr+(kc-1)*imr*jmr;
          ind3 = ic+alpha-1+(jc+beta-1)*imr+(kc-1)*imr*jmr;
          ind4 = ic-1+(jc-1)*imr+(kc+gamma-1)*imr*jmr;
          ind5 = ic+alpha-1+(jc-1)*imr+(kc+gamma-1)*imr*jmr;
          ind6 = ic-1+(jc+beta-1)*imr+(kc+gamma-1)*imr*jmr;
          ind7 = ic+alpha-1+(jc+beta-1)*imr+(kc+gamma-1)*imr*jmr;
          
          for (E_Int v = 4; v <= f.getNfld(); v++)
          {
            E_Float interp =
              cf[0] * f2(ind0,v)+
              cf[1] * f2(ind1,v)+
              cf[2] * f2(ind2,v)+
              cf[3] * f2(ind3,v)+
              cf[4] * f2(ind4,v)+
              cf[5] * f2(ind5,v)+
              cf[6] * f2(ind6,v)+
              cf[7] * f2(ind7,v);
            if (E_abs(f(ind,v)-interp) < E_abs((*error)(ind,v)))
              (*error)(ind,v) = f(ind,v)-interp;
          }
          if (f.getNfld() == 9 || f.getNfld() == 11)
            (*error)(ind,f.getNfld()) = f(ind,f.getNfld());
        }
        no++;
      }
    }
    errors.push_back(error);
  }
  cout << "done."<<endl;

  /*----------------------------------*/
  /* Checking non interpolated points */
  /*----------------------------------*/
  for (vector<FldArrayF*>::iterator it = errors.begin();
       it != errors.end(); it++)
  {
    FldArrayF& e = **it;
    
    for (E_Int ind = 0; ind < e.getSize(); ind++)
    {
      if ( e(ind,4) > 9.e5 ) // non interpolated point
      {
        for (E_Int v = 4; v <= e.getNfld(); v++)
        {
          e(ind,v) = 0.; // set to 0. and cellN to 0 also
        }
      }
    }
  }

  /*-----------------------*/
  /* Checking solid points */
  /*-----------------------*/
  E_Boolean chimera = E_False;
  if (errors[0]->getNfld() == 9 ||  errors[0]->getNfld() == 11)
    chimera = E_True;

  if (chimera == E_True)
  {
    E_Int cellN = errors[0]->getNfld();
    for (vector<FldArrayF*>::iterator it = errors.begin();
         it != errors.end(); it++)
    {
      FldArrayF& e = **it;
      
      for (E_Int ind = 0; ind < e.getSize(); ind++)
      {
        if ( e(ind,cellN) == 0 ) // solid point
        {
          for (E_Int v = 4; v < e.getNfld(); v++)
          {
            e(ind,v) = 0.; // set to 0.
          }
        }
      }
    }
  }
  
  /*----------*/
  /* L2 error */
  /*----------*/
  E_Int nfld = errors[0]->getNfld();
  FldArrayF L2error(nfld);
  L2error.setAllValuesAtNull();
  FldArrayF L0error(nfld);
  L0error.setAllValuesAtNull();

  for (E_Int v = 3; v < nfld; v++)
  {
    E_Int np = 0;
    for (vector<FldArrayF*>::iterator it = errors.begin();
         it != errors.end(); it++)
    {
      FldArrayF& e = **it;
      
      for (E_Int ind = 0; ind < e.getSize(); ind++)
      {
        L2error[v] = L2error[v] + e(ind,v+1)*e(ind,v+1);
        L0error[v] = E_max( L0error[v], e(ind,v+1) );
        np++;
      }
    }
    L2error[v] = sqrt(L2error[v]) / np;
  }
  
  E_Int errorsNb = 0;
  for (E_Int v = 3; v < nfld; v++)
  {
    if (L0error[v] >= 1.e-14)
    {
      cout << "INFO : L2 error ("<<v-2<<") = "<<L2error[v]<<endl;
      cout << "INFO : L0 error ("<<v-2<<") = "<<L0error[v]<<endl;
      errorsNb++;
    }
  }
  if (errorsNb == 0)
    cout << "INFO : Files are identical."<<endl;
  
  /*-----------------------*/
  /* Sauvegarde du fichier */
  /*-----------------------*/
  cout << "INFO : >> Saving error.plt..."<<flush;
  GenIO::getInstance()->tecwrite("error.plt", varString,
                                 imsm, jmsm, kmsm, errors);  
  cout << "done."<<endl;
}

//=============================================================================
// diff2
//=============================================================================
void diff2(E_Int argc, char** argv)
{
  char varString[256];
  char file1[256];
  char file2[256];
  strcpy(file1, argv[1]);
  strcpy(file2, argv[2]);
  
  /*---------------*/
  /* Reading file1 */
  /*---------------*/
  vector<FldArrayF*> field1;            // field read for each zone
  vector<E_Int> im1;
  vector<E_Int> jm1;
  vector<E_Int> km1;
  vector<FldArrayF*> ufield1;
  vector<FldArrayI*> c1;
  vector<E_Int> et1;
  cout << "INFO : << Reading "<<file1<<"..."<<flush;
  GenIO::getInstance()->tecread(file1, varString, field1, im1, jm1, km1, 
                                ufield1, c1, et1);
  cout << "done."<<endl;

  /*---------------*/
  /* Reading file2 */
  /*---------------*/
  vector<FldArrayF*> field2;            // field read for each zone
  vector<E_Int> im2;
  vector<E_Int> jm2;
  vector<E_Int> km2;
  vector<FldArrayF*> ufield2;
  vector<FldArrayI*> c2;
  vector<E_Int> et2;
  cout << "INFO : << Reading "<<file2<<"..."<<flush;
  GenIO::getInstance()->tecread(file2, varString, field2, im2, jm2, km2, 
                                ufield2, c2, et2);
  cout << "done."<<endl;

  /*-----------------------*/
  /* Computing error field */
  /*-----------------------*/
  vector<FldArrayF*> errors;
  vector<FldArrayF*>::iterator it2 = field2.begin();
  const E_Float EPS = 1.e-12;
  
  // We now try to identify the blocks by sizes and position
  for (vector<FldArrayF*>::iterator it = field1.begin();
       it != field1.end(); it++)
  {
    FldArrayF& f1 = **it;
    E_Int n1 = f1.getSize();
    FldArrayF* error = new FldArrayF(n1, f1.getNfld());
    error->setAllValuesAt(1.e6);
    errors.push_back(error);
    
    E_Boolean found = E_False;
    
    for (it2 = field2.begin();
         it2 != field2.end();
         it2++)
    {
      FldArrayF& f2 = **it2;
      E_Int n2 = f2.getSize();
      
      if ( n1 == n2 &&
           E_abs(f2(0,1) - f1(0,1)) < EPS &&
           E_abs(f2(0,2) - f1(0,2)) < EPS &&
           E_abs(f2(0,3) - f1(0,3)) < EPS &&
           E_abs(f2(n2-1,1) - f1(n1-1,1)) < EPS &&
           E_abs(f2(n2-1,2) - f1(n1-1,2)) < EPS &&
           E_abs(f2(n2-1,3) - f1(n1-1,3)) < EPS )
      {
        found = E_True;
        for (E_Int ind = 0; ind < n1; ind++)
        {
          E_Float x = f1(ind, 1);
          E_Float y = f1(ind, 2);
          E_Float z = f1(ind, 3);
          
          (*error)(ind, 1) = x;
          (*error)(ind, 2) = y;
          (*error)(ind, 3) = z;
          
          for (E_Int nf = 4; nf <= f1.getNfld(); nf++)
          {
            if (isnan(f1(ind, nf)))
            {
              cout << "FATAL : "<<file1<<" contains Nan values." <<endl;
              exit(0);
            }
            if (isnan(f2(ind, nf)))
            {
              cout << "FATAL : "<<file2<<" contains Nan values." <<endl;
              exit(0);
            }
            if (isinf(f1(ind, nf)))
            {
              cout << "FATAL : "<<file1<<" contains Infinite values." <<endl;
              exit(0);
            }
            if (isinf(f2(ind, nf)))
            {
              cout << "FATAL : "<<file2<<" contains Infinite values." <<endl;
              exit(0);
            }
            
            (*error)(ind, nf) = E_abs(f1(ind, nf) - f2(ind, nf));
          }
        }
        break;
      }
    } // iterator on field2
    if (found == E_False)
    {
      cout << "WARNING : A field on a block can not be compared."<<endl;
    } 
  }

  /*----------*/
  /* L2 error */
  /*----------*/
  E_Int nfld = errors[0]->getNfld();
  FldArrayF L2error(nfld);
  L2error.setAllValuesAtNull();
  FldArrayF L0error(nfld);
  L0error.setAllValuesAtNull();

  for (E_Int v = 3; v < nfld; v++)
  {
    E_Int np = 0;
    for (vector<FldArrayF*>::iterator it = errors.begin();
         it != errors.end(); it++)
    {
      FldArrayF& e = **it;
      
      for (E_Int ind = 0; ind < e.getSize(); ind++)
      {
        L2error[v] = L2error[v] + e(ind,v+1)*e(ind,v+1);
        L0error[v] = E_max( L0error[v], e(ind,v+1) );
        np++;
      }
    }
    L2error[v] = sqrt(L2error[v]) / np;
  }
  
  E_Int errorsNb = 0;
  for (E_Int v = 3; v < nfld; v++)
  {
    if (L0error[v] >= 1.e-14)
    {
      cout << "INFO : L2 error ("<< v-2 <<") = "<< L2error[v] <<endl;
      cout << "INFO : L0 error ("<< v-2 <<") = "<< L0error[v] <<endl;
      errorsNb++;
    }
  }
  if (errorsNb == 0)
    cout << "INFO : Files are identical."<<endl;

  /*-----------------------*/
  /* Sauvegarde du fichier */
  /*-----------------------*/
  cout << "INFO : >> Saving error.plt..."<<flush;
  GenIO::getInstance()->tecwrite("error.plt", varString, 
                                 im1, jm1, km1, errors);  
  cout << "done."<<endl;
  
}

//=============================================================================
// Main
//=============================================================================
int main(int argc, char** argv)
{
  if (argc < 3 || argc > 4)
  {
    cout << "CDiff : usage :"<<endl;
    cout << "---------------"<<endl;
    cout << "CDiff meshr.plt outputr.plt outputm.plt, where :"<<endl;
    cout << "meshr.plt : fine reference mesh (nodes)."<<endl;
    cout << "outputr.plt : fine reference solution (centers)."<<endl;
    cout << "outputm.plt : solution to be tested (centers)."<<endl;
    cout << "output file : error.plt"<<endl;
    cout << " -- or --" <<endl;
    cout << "CDiff output1.plt output2.plt, where :"<<endl;
    cout << "output1.plt and output2.plt are the same meshes."<<endl;
    exit(0);
  }
  
  cout << ">> Diffing files <<"<<endl;

  if (argc == 4)
  {
    // Diffing two solutions on two different meshes
    diff1(argc, argv);
  }
  else if (argc == 3)
  {
    // Diffing two solution on the same mesh
    diff2(argc, argv);
  }
}
