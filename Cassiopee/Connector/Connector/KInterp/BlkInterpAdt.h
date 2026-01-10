/*    
    Copyright 2013-2026 ONERA.

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

// ============================================================================
// @Name BlkInterpAdt
// @Memo Class of Alternative Digital Tree for search of interpolation cell.
// @See BlkInterpIcg
/* @Text

Validity
    This class manages an alternative digital tree on a curvilinear KMesh.
    This tree is then used by the method getInterpolationCell
    to find quickly the interpolation cell of a given point (the cell that
    contains the point). This class is an alternative to BlkInterpIcg.

Design

   There is no parameter for this class!
   Attention : en non structure tetra, les noeuds de l adt correspondent
   aux elements et non aux sommets des elts
*/
//=============================================================================
class BlkInterpAdt : public BlkInterpWithKMesh
{
  public: 

    ///+ 1- Constructor / Destructor
    
    /** 1- Destructor. */
    virtual ~BlkInterpAdt();
    void destroy();
    
    /** 2- Build the BlkInterpAdt on a BlkMesh.
        @param mesh       IN mesh on which is built the BlkInterpAdt.
    */
    BlkInterpAdt(KMesh& mesh);
    
    ///-
    
    ///+ 2- GET methods

    /* Recherche de la cellule d'interpolation pour un KMesh tetra
       IN: x,y,z: coordonnees du point a interpoler
       OUT: noelt: numero du tetra trouve
       OUT: cf: coeffs d interpolation tetra
       retourne 0 si pas trouve
    */
    virtual 
    short searchInterpolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int& noelt, FldArrayF& cf);

    /** 2- Get the index of interpolation cell and interpolation coefficients 
        for point (x,y,z) in the KMesh on which is built the BlkInterpAdt.
        The index of cell containing the point (x,y,z) (if found) is ic,jc,kc.
        Interpolation coefficients are stored in the cf[8] vector and should 
        be used to interpolate a function F defined on the cell edges by:
        F(x,y,z)=cf[0]*F(x0,y0,z0)+cf[1]*F(x1,y1,z1)+...
        where (x0,y0,z0)... are the coordinates of the cell with the 
        numerotation:
                        4   
                        o*************o 6
                     *  .          *  *
                  *     .       *     *
               *        .    *        *
          5 o*************o 7         *            K
            *           . *           *            .
            *           . *           *            .
            *           . *           *            .
            *         0 o.*...........o 2          ........ J
            *        .    *        *            .
            *     .       *     *            .
            *  .          *  *            I
          1 o*************o 3
    
        This function returns 1 if a cell containing the point (x,y,z) 
        exists, 0 otherwise.
    */
    virtual
    short searchInterpolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int &ic, E_Int &jc, E_Int &kc, 
                                  FldArrayF& cf);

    /** Same as before but for a vector of points. Vectorized version.
        @param coord    IN     coordinates of vector of (x,y,z) points to be interpolated
                                         only the part [ istart : iend ] is treated
        @param istart   IN     first computed point
        @param iend     IN     last computed point
        @param ic,jc,kc OUT    vector of index of cell containing point
        @param cf       OUT    vector of interpolation coefficients
        @param found    OUT    vector of boolean, true if point has been found
                               out vectors are modified only for [istart : iend]
    */
    virtual
    void searchInterpolationCellv(FldArrayF& coord, E_Int istart, E_Int iend,
                                  FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                                  FldArrayF& cf, FldArrayIS& found);
    
    virtual
    short searchInterpolationCellByJump(E_Float x, E_Float y, E_Float z,
                                        E_Int &ic, E_Int &jc, E_Int &kc, 
                                        FldArrayF& cf);
    virtual
    void searchInterpolationCellByJumpv(FldArrayF& coord,E_Int istart, E_Int iend,
                                        FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                                        FldArrayF& cf, FldArrayIS& found);
    
//     /** Same but for order 3 interpolations. */
//     virtual
//     short searchInterpolationCellHO( E_Float x, E_Float y, E_Float z,
//                                      E_Int& ic, E_Int& jc, E_Int& kc,
//                                      FldArrayF& cf );
//     virtual
//     void searchInterpolationCellHOv( FldArrayF& coord,
//                                      E_Int istart, E_Int iend,
//                                      FldArrayI& ic, FldArrayI& jc,FldArrayI& kc,
//                                      FldArrayF& cf, FldArrayIS& found );
   
    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) for structured meshes
        (nearest boundary cell) */
    virtual
    short getExtrapolationCell(E_Float x, E_Float y, E_Float z,
                               E_Int& ic, E_Int& jc, E_Int& kc,
                               FldArrayF& cf,
                               const FldArrayI& cellNatureField,
                               E_Int testNature,
                               E_Float& test,
                               InterpolationType interpType,
                               InterpMeshType interpMeshType=EXT_CENTERS,
                               E_Float cfMax=30.);

    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) on structured meshes */
    virtual short getExtrapolationCellStruct(E_Float x, E_Float y, E_Float z,
                                             FldArrayI& indi,
                                             FldArrayF& cf,
                                             E_Int order, E_Float cfMax,
                                             const FldArrayF& cellNatureField,
                                             InterpolationType interpType);
    
    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) for unstructured meshes */
    virtual 
    short getExtrapolationCellUnstr(E_Float x, E_Float y, E_Float z,
                                    E_Int& noelt,
                                    FldArrayF& cf,
                                    FldArrayI& indi,
                                    E_Int order,
                                    const FldArrayF& cellNatureField);

    virtual
    short getExtrapolationCoeffForCell(
      E_Float x, E_Float y, E_Float z,
      E_Int ic, E_Int jc, E_Int kc,
      FldArrayF& cf,
      const FldArrayI& cellNatureField,
      E_Int testNature, E_Int order,
      E_Int& is, E_Int& js, E_Int& ks, E_Float cfMax,
      InterpolationType interpType, 
      InterpMeshType interpMeshType=EXT_CENTERS);
    virtual
    short getExtrapolationCoeffForCell(
      E_Float x, E_Float y, E_Float z,
      E_Int ic, E_Int jc, E_Int kc,
      FldArrayF& cf, E_Float cfMax,
      InterpolationType interpType, 
      InterpMeshType interpMeshType=EXT_CENTERS);
    /* Recherche de la cellule d extrapolation pour un KMesh tetra
       IN: x,y,z: coordonnees du point a interpoler
       OUT: noelt: numero du tetra trouve
       OUT: cf: coeffs d interpolation tetra
       retourne 0 si pas trouve
    */
    virtual 
    short searchExtrapolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int& noelt, FldArrayF& cf, FldArrayI& ind,
                                  E_Int order, const FldArrayF& cellNatureField);
    virtual 
    short searchExtrapolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int& ic, E_Int& jc, E_Int& kc,
                                  FldArrayF& cf,
                                  E_Int order, E_Float cfMax,
                                  const FldArrayF& cellNatureField);

    ///-

  private:
    
    /* routine generale */
    void buildAdt();
    /* construit l'adt a partir d'un kmesh structure */
    void buildStructAdt();
    /* construit l'adt a partir d'un kmesh non structure */
    void buildUnstructAdt();
    
    void insert(E_Int ind, E_Float xmin, E_Float ymin, E_Float zmin,
                E_Float xmax, E_Float ymax, E_Float zmax);

    /* Recherche de la liste des cellules candidates. 
       Retourne la taille de listOfCandidateCells */
    E_Int getListOfCandidateCells(E_Float x, E_Float y, E_Float z,
                                  std::list<E_Int>& listOfCandidateCells);

    E_Float _xmax;   // bounding box of mesh
    E_Float _ymax;
    E_Float _zmax;
    E_Float _xmin;
    E_Float _ymin;
    E_Float _zmin;

    BlkIntTreeNode* _tree;
};

//=============================================================================

// ===== Interp/BlkInterpAdt.h === Last line ===
