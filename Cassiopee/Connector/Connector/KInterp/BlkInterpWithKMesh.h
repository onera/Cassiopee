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

// ============================================================================
/*
   This class manages interpolable blocks when the interpData is based on 
   a KMesh 
*/
// ============================================================================
class BlkInterpWithKMesh : public BlkInterpData
{
  public:
    
    ///+ 1- Constructor / Destructor
    /** Constructor */
    BlkInterpWithKMesh(KMesh& mesh);
   
    /** Destructor. */
    virtual
    ~BlkInterpWithKMesh();
    ///-
   
    /** Get extrapolation coefficients cf for point (x,y,z) inside cell 
        (ic,jc,kc) (in extended centers).
        
        If order is equal to zero, if at least an interpolable point 
        (cellNatureField=1) exists in this cell, return true and cf are 
        degenerated (order 0) extrapolation coefficients.
        
        If order is equal to one, try to find an interpolable tetrahedra. 
        Return true if OK, and cf are extrapolation cofficients in this 
        tetrahedra (order 1).

        If order is equal to two, we try first to find an interpolable 
        tetrahedra, if it is not possible we switch to an order 0 formula.
        
        CellNatureField is information relative to the grid we are 
        interpolating from.

        @param x,y,z               IN  coordinates of point to extrapolate
        @param ic,jc,kc            IN  index of cell that furnishes the values
        @param cf                  OUT extrapolation coefficients
        @param cellNatureField     IN  cell nature field on the grid we are extrapolating from (optional)
        @param order               IN  order of extrapolation formula (see up) (optional)
        @param is,js,ks            OUT if extrapolation or degeneration cannot be performed
                                       is,js,ks are the index of a neigbouring cell of ic,jc,kc
                                       nearer of (x,y,z) (optional).
     */
    virtual short getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                                               E_Int ic, E_Int jc, E_Int kc,
                                               FldArrayF& cf, E_Float cfMax,
                                               InterpolationType interpType,
                                               InterpMeshType interpMeshType=EXT_CENTERS);
    
    virtual short getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                                               E_Int ic, E_Int jc, E_Int kc,
                                               FldArrayF& cf,
                                               const FldArrayI& cellNatureField,
                                               E_Int testNature, E_Int order,
                                               E_Int& is, E_Int& js, E_Int& ks, E_Float cfMax,
                                               InterpolationType interpType,
                                               InterpMeshType interpMeshType=EXT_CENTERS);

    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) on structured meshes (nearest boundary cell) */
    virtual short getExtrapolationCell(E_Float x, E_Float y, E_Float z,
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

    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) on unstructured meshes */
    virtual short getExtrapolationCellUnstr(E_Float x, E_Float y, E_Float z,
                                            E_Int& noet,
                                            FldArrayF& cf,
                                            FldArrayI& indi,
                                            E_Int order,
                                            const FldArrayF& cellNatureField);

    /** Transformation de ic,jc,kc defini en centres etendus en centres
        Retourne 0 si extrapolation du bord, 1 sinon */
    virtual E_Int fromExtendedToStandardCenters( 
      E_Int ic, E_Int jc, E_Int kc,
      FldArrayI& indTab, InterpolationType interpType);
   
    /** Meme chose mais pour un ensemble de points entre istart et iend-1
        extrap : 0 si extrapolation du bord, 1 si interieur */
    virtual void fromExtendedToStandardCentersv( 
      E_Int istart, E_Int iend, 
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayI& indTab, FldArrayI& extrap, InterpolationType interpType);  
  
    /* Calcul des indices i,i+1... de la molecule d' interpolation 
       Le premier element de indi : indi[0] donne l info sur l ordre 
       d interpolation et n est pas calcule ici*/
    virtual void compStandardIndices(
      E_Int ic, E_Int jc, E_Int kc, 
      FldArrayI& indi, InterpolationType interpType);
    /* Calcul des indices i,i+1... des molecules d' interpolation 
       comprises entre istart et iend-1 inclus
       Les premiers elemts de indiv :indiv(i,1) sont deja calcules avant
       et valent l ordre d interpolation local au pt  */
    virtual void compStandardIndicesv(
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc, 
      FldArrayI& indi, InterpolationType interpType);
    ///-
  protected:
    /** Interpolation d ordre 2 non structuree */
    virtual 
    short searchInterpolationCellUnstruct(E_Float x, E_Float y, E_Float z,
                                          E_Int& noelt,
                                          FldArrayI& indi, FldArrayF& cf);

    /** 2nd order interpolation search using O2CF interpolation type : 8 coefs. */
    virtual
    short searchInterpolationCellO2CF( E_Float x, E_Float y, E_Float z,
                                       E_Int& ic, E_Int& jc, E_Int& kc,
                                       FldArrayF& cf ) ;
    /** 2nd order interpolation search using O2ABC interpolation type. */
    virtual 
    short searchInterpolationCellO2ABC(E_Float x, E_Float y, E_Float z, 
                                       E_Int& ic, E_Int& jc, E_Int& kc);
    
    /** 3rd order interpolation not relevant for Ext Centers Interp.
     Returns -1*/
    virtual
    short searchInterpolationCellO3CF( E_Float x, E_Float y, E_Float z,
                                       E_Int& ic, E_Int& jc, E_Int& kc,
                                       FldArrayF& cf );

    /** 3rd order interpolation search using O3ABC interpolation type. */
    virtual short searchInterpolationCellO3ABC(
      E_Float x, E_Float y, E_Float z, E_Int& ic, E_Int& jc, E_Int& kc);

    /** 5th order interpolation search using O5ABC interpolation type. */
    virtual short searchInterpolationCellO5ABC(
      E_Float x, E_Float y, E_Float z, E_Int& ic, E_Int& jc, E_Int& kc);
    
    /** O2CF interpolation search. Vectorized. */
    virtual void searchInterpolationCellO2CFv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayF& cf, FldArrayIS& found);
   
    /** O2ABC interpolation search. Vectorized. */
    virtual void searchInterpolationCellO2ABCv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayIS& found);
   
    /** O3CF interpolation search. Vectorized.
        O3CF is only done for InterpCart type. Exits for all the others 
        interpolation data types (with extended centers).*/
    virtual void searchInterpolationCellO3CFv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayF& cf, FldArrayIS& found);

    /** O3ABC interpolation search. Vectorized. */
    virtual void searchInterpolationCellO3ABCv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayIS& found);

    /** O5ABC interpolation search. Vector. */
    virtual void searchInterpolationCellO5ABCv( 
      FldArrayF& coord, E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayIS& found);

    /** Compute the Lagrange coefs for interpolation of (x,y,z)
        Retourne 1 si l ordre est rabaisse a 2 au pt */
    virtual short compLagrangeCoefs(
      E_Float x, E_Float y, E_Float z,
      E_Int ic, E_Int jc, E_Int kc,
      FldArrayF& cf,
      InterpolationType interpType,
      InterpMeshType interpMeshType=EXT_CENTERS);

    virtual
    short searchInterpolationCellByJump(E_Float x, E_Float y, E_Float z,
                                        E_Int& ic, E_Int& jc, E_Int& kc,
                                        FldArrayF& cf);

    virtual 
    void searchInterpolationCellByJumpv( FldArrayF& coord,
                                         E_Int istart, E_Int iend,
                                         FldArrayI& ic, FldArrayI& jc,FldArrayI& kc,
                                         FldArrayF& cf, FldArrayIS& found );
    
    // methodes purement virtuelles 
    // non structure
    virtual 
    short searchInterpolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int& noelt, FldArrayF& cf) = 0;
    // structure
    virtual
    short searchInterpolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int &ic, E_Int &jc, E_Int &kc, 
                                  FldArrayF& cf) = 0;
    virtual
    void searchInterpolationCellv(FldArrayF& coord, E_Int istart, E_Int iend,
                                  FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
                                  FldArrayF& cf, FldArrayIS& found) = 0;
    virtual 
    short searchExtrapolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int& noelt, FldArrayF& cf, FldArrayI& indi, 
                                  E_Int order, const FldArrayF& cellNatureField) = 0;
    virtual 
    short searchExtrapolationCell(E_Float x, E_Float y, E_Float z,
                                  E_Int& ic, E_Int& jc, E_Int& kc,
                                  FldArrayF& cf,
                                  E_Int order, E_Float cfMax,
                                  const FldArrayF& cellNatureField) = 0;

    void coordHexa(E_Int ind, E_Int ni, E_Int nj,
                   E_Float* xl, E_Float* yl, E_Float* zl,
                   E_Int &ic, E_Int &jc, E_Int &kc,
                   E_Float* xt, E_Float* yt, E_Float* zt);
    void coordHexal(E_Int ind, E_Int ni, E_Int nj,
                    E_Float* xl, E_Float* yl, E_Float* zl,
                    E_Int &ic, E_Int &jc, E_Int &kc,
                    E_Float* xt, E_Float* yt, E_Float* zt);
    
    KMesh& _mesh;      // Mesh on which is built interpWithKMesh
};

// ======== Connector/Interp/BlkInterpWithKMesh.h ============================
