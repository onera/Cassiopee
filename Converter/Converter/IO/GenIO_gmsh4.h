  if (nBAR > 0) {
    cnBAR = new FldArrayI(nBAR, 2, true, fo);
    connect.push_back(cnBAR);
    eltType.push_back(1); nBAR = 0; }  
  if (nTRI > 0) {
    cnTRI = new FldArrayI(nTRI, 3, true, fo);
    connect.push_back(cnTRI);
    eltType.push_back(2); nTRI = 0; }
  if (nQUAD > 0) {
    cnQUAD = new FldArrayI(nQUAD, 4, true, fo);
    connect.push_back(cnQUAD);
    eltType.push_back(3); nQUAD = 0; }
  if (nTETRA > 0) {
    cnTETRA = new FldArrayI(nTETRA, 4, true, fo);
    connect.push_back(cnTETRA);
    eltType.push_back(4); nTETRA = 0; }
  if (nHEXA > 0) {
    cnHEXA = new FldArrayI(nHEXA, 8, true, fo);
    connect.push_back(cnHEXA);
    eltType.push_back(7); nHEXA = 0; }
  if (nPENTA > 0) {
    cnPENTA = new FldArrayI(nPENTA, 6, true, fo);
    connect.push_back(cnPENTA);
    eltType.push_back(6); nPENTA = 0; }
  if (nPYRA > 0) {
    cnPYRA = new FldArrayI(nPYRA, 5, true, fo);
    connect.push_back(cnPYRA);
    eltType.push_back(5); nPYRA = 0; }
  if (nNODE > 0) {
    indNODE = new FldArrayI(nNODE);
    connect.push_back(new FldArrayI());
    eltType.push_back(0);
    nNODE = 0; }
  if (nBAR_3 > 0) {
    cnBAR_3 = new FldArrayI(nBAR_3, 2, true, fo);
    connect.push_back(cnBAR_3);
    eltType.push_back(1); nBAR_3 = 0; }
  if (nTRI_6 > 0) {
    cnBAR_3 = new FldArrayI(nTRI_6, 2, true, fo);
    connect.push_back(cnTRI_6);
    eltType.push_back(1); nTRI_6 = 0; }
  if (nQUAD_8 > 0) {
    cnQUAD_8 = new FldArrayI(nQUAD_8, 2, true, fo);
    connect.push_back(cnQUAD_8);
    eltType.push_back(1); nQUAD_8 = 0; }
  if (nQUAD_9 > 0) {
    cnQUAD_9 = new FldArrayI(nQUAD_9, 2, true, fo);
    connect.push_back(cnQUAD_9);
    eltType.push_back(1); nQUAD_9 = 0; }
  if (nTETRA_10 > 0) {
    cnTETRA_10 = new FldArrayI(nTETRA_10, 2, true, fo);
    connect.push_back(cnTETRA_10);
    eltType.push_back(1); nTETRA_10 = 0; }
  if (nPYRA_14 > 0) {
    cnPYRA_14 = new FldArrayI(nPYRA_14, 2, true, fo);
    connect.push_back(cnPYRA_14);
    eltType.push_back(1); nPYRA_14 = 0; }
  if (nPENTA_15 > 0) {
    cnPENTA_15 = new FldArrayI(nPENTA_15, 2, true, fo);
    connect.push_back(cnPENTA_15);
    eltType.push_back(1); nPENTA_15 = 0; }
  if (nHEXA_20 > 0) {
    cnHEXA_20 = new FldArrayI(nHEXA_20, 2, true, fo);
    connect.push_back(cnHEXA_20);
    eltType.push_back(1); nHEXA_20 = 0; }
  if (nHEXA_27 > 0) {
    cnHEXA_27 = new FldArrayI(nHEXA_27, 2, true, fo);
    connect.push_back(cnHEXA_27);
    eltType.push_back(1); nHEXA_27 = 0; }
