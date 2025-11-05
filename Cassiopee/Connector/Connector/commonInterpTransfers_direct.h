imdjmd           = imd*jmd;
E_Int max_thread = min(nvars , E_Int(__NUMTHREADS__));


#pragma omp parallel default(shared) num_threads(max_thread)
 {

#ifdef _OPENMP
  E_Int  ithread           = omp_get_thread_num()+1;
  E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
  E_Int ithread = 1;
  E_Int Nbre_thread_actif = 1;
#endif
  // Calcul du nombre de champs a traiter par chaque thread
  E_Int chunk = nvars/Nbre_thread_actif;
  E_Int r = nvars - chunk*Nbre_thread_actif;
  E_Int eq_deb, eq_fin;
  // equations traitees par thread
  if (ithread <= r) 
  { eq_deb = (ithread-1)*(chunk+1); eq_fin = eq_deb + (chunk+1); }  
  else { eq_deb = (chunk+1)*r+(ithread-r-1)*chunk; eq_fin = eq_deb + chunk; }

  E_Float* ptrCoefs = donorCoefsF->begin();
  E_Int indR, type;
  E_Int indD0, indD, ncfLoc;
  E_Int noi = 0; // compteur sur le tableau d indices donneur
  E_Int sizecoefs = 0;

  for (E_Int noind = 0; noind < nbRcvPts; noind++)
  { 
    //
    // adressage direct pour indR
    //
    indR = noind;
# include "commonInterpTransfers.h" 
    ptrCoefs += sizecoefs;
  } 

 }// omp
