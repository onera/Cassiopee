//=============================================================================
// This routine reads an output_skin file 
// and convert it to the generic USURP format
//=============================================================================

#define _ISO_VECTOR_
# include "usurp/Convert/main.h"

using namespace E_FUNC;

extern "C"
{
  void writebcfile_(const E_Int& nzone, const E_Int* nit, const E_Int* njt, 
                    const E_Int* nkt);
  void writeibfile_(const E_Int& nzone, const E_Int& ncellmax, 
                    const E_Int* ncelltot, const E_Int* iblank);
  
  void writegrid_(const E_Int& nzone, const E_Int& nptsmax,
                  const E_Int* nit,const E_Int* njt, const E_Int* nkt,
                  const E_Float* xb, const E_Float* yb, const E_Float* zb);
}

//=============================================================================
/* Convert the cellnaturefield defined at nodes here in centers
   Respect the OVERFLOW syntax
   Cell nature field is supposed to be last field.
*/
//=============================================================================
void writeIBlankArray(FldArrayI& nit, FldArrayI& njt, 
                      vector<FldArrayF*>& field)
{
  E_Int nzone = field.size();
  E_Int nfld = (field[0])->getNfld();
  E_Int ncellmax = 0;
  E_Int ind1, ind2, ind3, ind4;
  E_Int indCell;
  E_Int prod;
  E_Int ni, ni1,nj1, ni1nj1;
  E_Float iv1, iv2, iv3, iv4;
  FldArrayI ncelltot(nzone);

  for ( E_Int n = 0; n < nzone; n++)
  { 
    ni1nj1 = (nit[n]-1)*(njt[n]-1);
    ncellmax = E_max(ncellmax, ni1nj1);
    ncelltot[n] = ni1nj1;
  }

  FldArrayI iblank(ncellmax, nzone);
  iblank.setAllValuesAt(-10);

  for(E_Int n = 0; n < nzone; n++)
  {
    ni = nit[n];
    ni1 = ni-1;
    nj1 = njt[n]-1;
    FldArrayF& oneField = *(field[n]);
    for (E_Int j = 1; j <= nj1; j++)
      for (E_Int i = 1; i <= ni1; i++)
      {
        ind1 = (i-1) + (j-1)*ni;
        ind2 = ind1 + 1;
        ind3 = ind2 + ni;
        ind4 = ind1 + ni;
        indCell = (i-1) + (j-1)*ni1;
        
        iv1 = oneField(ind1, nfld);
        iv2 = oneField(ind2, nfld);
        iv3 = oneField(ind3, nfld);
        iv4 = oneField(ind4, nfld);
        
        prod = E_Int(iv1*iv2*iv3*iv4);
        if ( prod == 0)
          iblank(indCell, n+1) = 0;
        else if (prod == 8)
          iblank(indCell, n+1) = -1;
        else 
          iblank(indCell, n+1) = 1;
      }
  }

  writeibfile_(nzone, ncellmax, ncelltot.begin(), iblank.begin());

}

//=============================================================================
// main
//=============================================================================
int main(E_Int argc, char** argv)
{
  char infile[256];
  char varString[256];
  if (argc <= 1 || argc > 2)
  {
    cout << "tec2usurp : convert binary tecplot skin file to USURP format "<<endl;
    cout << "usage : tec2usurp <output_skin.plt> "<<endl;
    cout << "<output_skin.plt> : overlapped surfaces (in)"<<endl;
    cout << "output : generic.ib, generic.grd, generic.bc."<<endl;
    exit(0);
  }
  
  strcpy(infile, argv[1]);
  
  /*----------------------------------------------------------*/
  /* Reading output_skin.plt and store fields in structBlocks */
  /*----------------------------------------------------------*/
  vector<E_Int> ni;
  vector<E_Int> nj;
  vector<E_Int> nk;
  vector<FldArrayF*> f;
  vector<FldArrayI*> c;
  vector<FldArrayF*> uf;
  vector<E_Int> et;

  cout << "Reading..."<<flush;
  if (GenIO::getInstance()->tecread(infile, varString, f, ni, nj, nk,
                                    uf, c, et) != 0)
  {
    cout << "FATAL : tec2usurp : can not open file"<<endl; exit(0);
  }
    
  E_Int nzone = f.size();

  cout << " ("<<nzone<<" zones read)."<<endl;
  
  FldArrayI nit(nzone);
  FldArrayI njt(nzone);
  FldArrayI nkt(nzone);

  for (E_Int n = 0; n < nzone; n++)
  {
    nit[n] = ni[n];
    njt[n] = nj[n];
    nkt[n] = nk[n];
  }
 
  // Write the generic.bc file with data for meshes
  cout << "Writing generic.bc..."<<flush;
  writebcfile_(nzone, nit.begin(), njt.begin(), nkt.begin());
  cout << "done."<<endl;

  // Write iblank file for all grids
  // first convert in centers the cellN to iblank
  cout << "Writing generic.ib..."<<flush;
  writeIBlankArray(nit, njt, f);
  cout <<"done."<<endl;

  //----------------------------------------------
  // Write the grid file generic.grd
  //----------------------------------------------
  E_Int nptsmax = 0;
  for (E_Int n = 0; n < nzone; n++)
    nptsmax = E_max(nptsmax, ni[n]*nj[n]*nk[n]);
   
  FldArrayF coordx(nptsmax, nzone); //max allocated for fortran
  FldArrayF coordy(nptsmax, nzone);
  FldArrayF coordz(nptsmax, nzone);
  for (E_Int n = 0; n < nzone; n++)
  {
    FldArrayF& field = *f[n];
    for(E_Int i = 0; i < field.getSize(); i++)
    {
      coordx(i, n+1) = field(i, 1);
      coordy(i, n+1) = field(i, 2);
      coordz(i, n+1) = field(i, 3);
    }
  }
  cout << "Writing generic.grd..."<<flush;
  writegrid_(nzone, nptsmax, nit.begin(), njt.begin(), nkt.begin(),  
             coordx.begin(), coordy.begin(), coordz.begin());
  cout << "done."<<endl;

  ni.clear();
  nj.clear();
  nk.clear();
  f.clear();

  return 1;
}




