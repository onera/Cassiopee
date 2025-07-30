  /*------------------------------------------*/
  /* Extraction des coordonnees du pt corrige */
  /*------------------------------------------*/
  FldArrayF* coordxPC; FldArrayF* coordyPC; FldArrayF* coordzPC;
  E_Int okc1 = K_NUMPY::getFromNumpyArray(pyArrayXPC, coordxPC);
  E_Int okc2 = K_NUMPY::getFromNumpyArray(pyArrayYPC, coordyPC);
  E_Int okc3 = K_NUMPY::getFromNumpyArray(pyArrayZPC, coordzPC);
  E_Float* xPC = coordxPC->begin();
  E_Float* yPC = coordyPC->begin();
  E_Float* zPC = coordzPC->begin();

  /*----------------------------------------*/
  /* Extraction des coordonnees du pt paroi */
  /*----------------------------------------*/
  FldArrayF* coordxPW;
  FldArrayF* coordyPW;
  FldArrayF* coordzPW;
  E_Int okw1 = K_NUMPY::getFromNumpyArray(pyArrayXPW, coordxPW);
  E_Int okw2 = K_NUMPY::getFromNumpyArray(pyArrayYPW, coordyPW);
  E_Int okw3 = K_NUMPY::getFromNumpyArray(pyArrayZPW, coordzPW);
  E_Float* xPW = coordxPW->begin();
  E_Float* yPW = coordyPW->begin();
  E_Float* zPW = coordzPW->begin();

  /*--------------------------------------------*/
  /* Extraction des coordonnees du pt interpole */
  /*--------------------------------------------*/
  FldArrayF* coordxPI; FldArrayF* coordyPI; FldArrayF* coordzPI;
  E_Int oki1 = K_NUMPY::getFromNumpyArray(pyArrayXPI, coordxPI);
  E_Int oki2 = K_NUMPY::getFromNumpyArray(pyArrayYPI, coordyPI);
  E_Int oki3 = K_NUMPY::getFromNumpyArray(pyArrayZPI, coordzPI);
  E_Float* xPI = coordxPI->begin();
  E_Float* yPI = coordyPI->begin();
  E_Float* zPI = coordzPI->begin();

  /*--------------------------------------------*/
  /* Extraction des var post traitement */
  /*--------------------------------------------*/
  FldArrayF* densF; 
  // FldArrayF* pressF; FldArrayF* utauF; FldArrayF* yplusF;
  // FldArrayF* vxF; FldArrayF* vyF; FldArrayF* vzF; FldArrayF* kcurvF;

  E_Int okD = K_NUMPY::getFromNumpyArray(pyArrayDens, densF);

  E_Float* density  = densF->begin();
  // E_Float* pressure = densF->begin()+nbRcvPts;
  // E_Float* vx = densF->begin()+nbRcvPts*2;
  // E_Float* vy = densF->begin()+nbRcvPts*3;
  // E_Float* vz = densF->begin()+nbRcvPts*4;
