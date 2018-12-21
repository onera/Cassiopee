  if (nBAR > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nTRI > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nQUAD > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nTETRA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nHEXA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nPENTA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nPYRA > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nNODE > 0)
  {
    FldArrayF* an = new FldArrayF(nNODE,3);
    for (E_Int i = 0; i < nNODE; i++) 
    { ind = (*indNODE)[i]-1; ind = indirNodes[ind]-1;
      (*an)(i,1) = f1[ind]; (*an)(i,2) = f2[ind]; (*an)(i,3) = f3[ind]; }
    unstructField.push_back(an);
    delete indNODE;
  }
  if (nTRI_6 > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
  if (nQUAD_9 > 0)
  {
    FldArrayF* an = new FldArrayF(f);
    unstructField.push_back(an);
  }
