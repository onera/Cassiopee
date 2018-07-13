  /*-------------------------------------*/
  /* Extraction des indices des donneurs */
  /*-------------------------------------*/
  FldArrayI* donorPtsI;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyIndDonor, donorPtsI, true);
  E_Int* donorPts = donorPtsI->begin();
  /*----------------------*/
  /* Extraction des types */
  /*----------------------*/
  FldArrayI* typesI;
  E_Int res_type=  K_NUMPY::getFromNumpyArray(pyArrayTypes, typesI, true);
  E_Int* types  = typesI->begin();
  E_Int nbRcvPts= typesI->getSize();// taille du numpy = nb de pts a interpoler
  /*-----------------------*/
  /* Extraction des coefs  */
  /*-----------------------*/
  FldArrayF* donorCoefsF;
  E_Int res_coef = K_NUMPY::getFromNumpyArray(pyArrayCoefs, donorCoefsF, true);
