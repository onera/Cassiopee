/*    
    Copyright 2013-2019 Onera.

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
# ifndef _POST_ZIPPER_STRUCT_BLOCK_H_
# define _POST_ZIPPER_STRUCT_BLOCK_H_

# include "kcore.h"
# include <vector>
# include <list>
# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI
# define FldArrayIS K_FLD::FldArrayIS
# define FldArrayB K_FLD::FldArrayB
//=============================================================================
/* Class defining a structured block for zipper */
//=============================================================================
class StructBlock
{
  public: 
    ///+ 1- Constructors / Destructor
    
    /** Constructor. Define a structured block.
        IN: id: block number,
        IN: ni,nj,nk: the number of points in the 3 directions,
        IN: field: field made of coord + fields. 
        IN: posx, posy, posz: position of coordinates in field
        IN: overlapTol: overlap tolerance 
        IN: matchTol: matching bcs tolerance
    */
    StructBlock(E_Int id, E_Int ni, E_Int nj, E_Int nk,
                E_Int posx, E_Int posy, E_Int posz,
                E_Float overlapTol, E_Float matchTol,
                FldArrayF& field);

    /** Destructor */
    ~StructBlock();

    /** Get the match tolerance */
    E_Float getMatchTol();

    /** Get mesh data, number of points, coordinates */
    E_Int getIm();
    E_Int getJm();
    E_Int getKm();
    FldArrayF& getCoord();

    /** Get the cfdField attached to StructBlock */
    FldArrayF& getCfdField();

    /** Get the global field ( (x,y,z) + cfdField ) attached to StructBlock */
    FldArrayF& getGlobalField();

    /** Get the block Id */
    E_Int getBlockId();

    /** Get the number of mesh points */
    E_Int getNbOfMeshPts(); 

    /** Get the number of CFD fields */
    E_Int getNbOfCfdFields();
    
    /** Get the number of fields (x,y,z)+CFD fields */
    E_Int getNbOfFields();

    /** Return the Iblank Array */
    FldArrayIS&  getIBlankArray();
    
    /** Return the strings of block */
    std::list<FldArrayI*>& getStrings();

    /** Return the points to be stringed */
     FldArrayI& getPointsForStrings();

    /** Compute the boundary nodes of the StructBlock, the degenerations,
        and the triangular connectivity */
    void compIBndIDgUnsConnectEN();

    /** Return the list of nodes to be merged */
    FldArrayI& getIBnd();
    
    /** Return the size of the list of nodes to be merged */
    E_Int getIBndSize();

    /** Return the array of degenerations of nodes */
    FldArrayI& getIDg();

    /** Return the size of the array of degenerations of nodes */
    E_Int getIDgSize();

    /** Return the array of triangular connectivity */
    FldArrayI& getUnsConnectEN();

    /** Return the size of the array of triangular connectivity */
    E_Int getUnsConnectENSize();

    /** Increment the indices of nodes in _ibnd in order to copy _ibnd
        in global ibndG */
    void incrIBnd(E_Int incr);

    /** Increment the indices of nodes in _idg in order to copy _ibnd
        in global idgG */
    void incrIDg(E_Int incr);
    
    /** Increment the indices of nodes in _unsConnectEN in order to copy
        _unsConnectEN in global unsConnectENG */
    void incrUnsConnectEN(E_Int incr);
    
    /** Compute iblank field for a given block 
        optionally: 
        dir specifies a blanking in a direction: 1:i, 2:j, 3:k
        from index imin to imax */
    void compIBlank(std::vector<StructBlock*>& vectOfBlks, E_Int noBlk2);
                  
    /** Compute iblank field for cell centers  for a given block */
    void compIBlankCC();
    
    /** Returns the indirection array */
    FldArrayI& getIndirectionArray();
    
    /** Returns the interpolable nodes array */
    FldArrayI& getInterpolableNodeArray();

    /** Returns the interpolation nodes array */
    FldArrayI& getInterpolationNodeArray();

    /** Return the number of interpolable pts with pts of block noBlk */
    E_Int getNbOfInterpolablePtsForBlk(E_Int noBlk);
    
    /** Computes information for interpolable cells */
    void selectTypeOfBlks(std::vector<StructBlock*>& vectOfBlks);
    
    /** Return the list of overlapping blocks */
    std::vector<StructBlock*>& getListOfOverlappingBlks();

    /** Return the list of matching blocks */
    std::vector<StructBlock*>& getListOfMatchingBlks();
    
    /** Clean lists of blocks */
    void cleanListsOfBlks( std::vector<StructBlock*>& vectOfBlks);
    
    /** Update the list of overlapping blocks with respect to blanking */
    void addOverlappingBlk(StructBlock* overlappingBlk);
    
    /** Check Interpolation information. If there is a point closer to the
       interpolable pts of current block that come from a nonoverlapping blk,
       erase the interpolation blk of the list _overlappingBlk.
       Useful for NS matching cells */
    void checkValidityOfInterpolationCell(
      E_Int noBlk,
      std::vector<StructBlock*>& vectOfBlks);
    
    /** Computes information for interpolable cells */
    void compInterpolationInformation(std::vector<StructBlock*>& vectOfBlks);
    
    /** Write the block in unstructured format, considering iblank.
     add = true, appends the block. */
    void write(char* fileName, E_Boolean add = false);

    /** Write indices of points in line shape */
    void writeLine(FldArrayI& indices,
                   char* fileName, E_Boolean add = false);
    
    /** Make 1D strings of indices of pts identified for stringing */
    void stringing(std::list<FldArrayI*>& strings);

    /** Find the points used in the stringing.
     noBlk : number of current block
     vectOfBlocks : all blocks vector. */
    void identifyGapBndPts( E_Int noBlk,
                            std::vector<StructBlock*> vectOfBlocks);
    
    /** Update the iblank array for points that are in an overlap but were
        not blanked at the 1st time */
    void updateIBlankArray();

    /* Compute the links between points.
       inIndices is an array of indices of mesh points.
       links is an array telling if a point is linked to another
       point of the inIndices array. Points are said to be linked if
       they are neighbour on the grid */
    void compLinks(FldArrayI& inIndices, FldArrayI& links);
   
    /* Chain points in inIndices */
    E_Int chainPoints(FldArrayI& inIndices, FldArrayI& links, 
                      FldArrayI& string, FldArrayB& dejaVu);
    
  private:
 /* Projection du point (x,y,z) sur la frontière dir (i=1...). 
    Retourne false si le projeté n'est pas situé sur cette frontiere, 
    a un epsilon près.*/
    E_Boolean projectOrtho(E_Float x, E_Float y, E_Float z, E_Int dir);
                             
    /* Given a pt (x1,y1,z1), says if it belongs to a match/nearmatch/nomatch 
       bnd (WITH cellNF = 1 criterion for matching bnd)*/
    E_Boolean searchForMatchingBnd(E_Float x1, E_Float y1, E_Float z1,
                                   StructBlock* blk2);
    /* Given a pt (x1,y1,z1) says if it has a matching pt on blk2 
      WITHOUT cellNF = 1 criterion*/
    E_Boolean searchForMatchingBnd2(E_Float x1, E_Float y1, E_Float z1,
                                    StructBlock* blk2, E_Int& ind2);
    
    /* Test the intersection of bounding boxes of block1 and block2 */
    E_Boolean testBBIntersection(StructBlock& block1, StructBlock& block2);

    /* Return the bounding box of mesh */
    void getBoundingBox(E_Float& xmax, E_Float& ymax, E_Float& zmax,
                        E_Float& xmin, E_Float& ymin, E_Float& zmin);
    
    /* Return the bounding box of cell i of mesh */
    void getBoundingBoxOfCell(E_Int i,
                              E_Float& xmax, E_Float& ymax, E_Float& zmax,
                              E_Float& xmin, E_Float& ymin, E_Float& zmin);

    /* Compute the distance of a node of index ind to the cell of first index
       indCell2. Return the minimum distance between ind and each vertex of
       the cell of index indCell2*/
    void computeDistanceOfNodeToCell(std::vector<StructBlock*>& vectOfBlks,
                                     E_Int ind,
                                     E_Float x, E_Float y, E_Float z,
                                     E_Int noBlk2, E_Int indCell2,
                                     E_Float& dist);
    
    /* Test if two blocks have a matching border */
    E_Boolean testMatchingBlks(StructBlock* blk1, StructBlock* blk2);
    
    /* Eliminate indices that are defined more than once in incidices
       array */
    void eliminateDoublePoints(FldArrayI& indices);

    /* Computes the interior points to be stringed */
    void compInteriorPtsForStringing(E_Int im, E_Int jm, E_Int& n);
    
    /* Compute the boundary points to be stringed */
    void compBndPtsForStringing(E_Int im, E_Int jm, E_Int& n);

    /* Compute the corners to be stringed */
    void compCornersForStringing(E_Int im, E_Int jm, E_Int& n);

    /* Identify the overlapping borders for stringing */
    void compOverlapBordersForString( E_Int& n );
    
    /* Identify the matching borders for stringing */
    void compMatchingBordersForString( E_Int noBlk, E_Int& n,
                                       std::vector<StructBlock*>& vectOfBlks );
    
    /* Test the Validity of the link between ind1 and ind2 in the direction
       dir. If they are connected by a line that is inside the mesh resulting
       of the blanking then return false*/
    E_Boolean testValidityOfLink(E_Int i, E_Int j, E_Int k,
                                 E_Int dir, E_Int sens);
    
    /* Compute Links for Interior Pts */
    void compLinksForInteriorPts(FldArrayI& inIndices, E_Int ir, 
                                 E_Int i, E_Int j, E_Int k,
                                 FldArrayI& links);
    
    /* Compute Links for Interior Pts */
    void compLinksForBndPts(FldArrayI& inIndices, E_Int ir, 
                            E_Int i, E_Int j, E_Int k,
                            FldArrayI& links);
    /* Test if the border node of index ind has to be set in a string
       or not */
    void testForBlankedBorderNode(E_Int ind, E_Int indCell1,
                                  E_Int indCell2,
                                  E_Int& n);

    /* Compute the list of blocks that have a matching bnd with a noBlk blk 's
       bnd */
    std::vector<StructBlock*>
    compListOfMatchingBlks(std::vector<StructBlock*>& vectOfBlks);  
    
    /* Given a vector of blks, eliminate double elements  */
    void eliminateDoubleElts(std::vector<StructBlock*>& vectOfBlks);
    /* Given a pt (x1,y1,z1) says if it has a matching pt on blk2 */
    E_Boolean isAMatchingBnd(E_Float x1, E_Float y1, E_Float z1,
                             StructBlock* blk2);

    /* Given a direction dir, return the array of indices of pts on the
       corresponding bnd
       dir = 1 : border i = 1  
             2 : border i = im
             3 : border j = 1
             4 : border j = jm */
    void createArrayOfBorders(E_Int dir, FldArrayI& indBord);

    /* Given a segment [ind11, ind12] on a boundary of direction dir2 
       and a matching block blk2, add the segment [ind11, ind12] in the string
       list. Return false if the segment is not added in the string. 
       out : n : index of next element of the string.*/
    E_Boolean addPtsInString(E_Int ind11, E_Int ind12, E_Int dir2,
                             StructBlock* blk2, E_Int& n);
    
    /* Compute the boundary nodes of the StructBlock, taking into account
       blanking (IBlank) */
    void compIBnd();

    /* Seek for boundary points among interior nodes */
    void compIBndInt(E_Int im, E_Int jm, E_Int& n);

    /* Seek for boundary points among boundary nodes */
    void compIBndExt(E_Int im, E_Int jm, E_Int& n);

    /* Compute the degenerations */
    void compIDg();
    
    /* Update _ibnd taking into account the degenerations */
    void updateIBnd();

    /* Compute the triangular connectivity taking into account
       the degenerations and the boundaries. */
    void compUnsConnectEN();

    /* Compute bounding box of cell */
    void boundingBoxOfCell(
      E_Int ind,
      E_Float& xmax, E_Float& ymax, E_Float& zmax, 
      E_Float& xmin, E_Float& ymin, E_Float& zmin);

    /* Return true if string is a loop (starting point = end point) */
    E_Boolean isStringALoop(FldArrayI& string);

  private:
    /* Structured mesh of nodes */
    E_Int _im;
    E_Int _jm;
    E_Int _km;
    FldArrayF _coord;

    /* CFD field on nodes */
    FldArrayF _cfdField;
    FldArrayF _globalField;

    /* Number of Cfd fields */
    E_Int _nCfdFields;

    /* nfld: nb of fields for cfdField array */
    E_Int _nfld;
    
    /* Number of fields ( (x,y,z) + Cfd fields ) */
    E_Int _nFields;

    /* Block Id */
    E_Int _id;
    
    /* iblank : is cell vertex in an overlap? */
    FldArrayIS _iblank;

    /* iblankCC : is cell center in an overlap? */
    FldArrayIS _iblankCC;

    /* Indirection tab */
    FldArrayI _indirection;

    /* Interpolable nodes of current block */
    FldArrayI _interpolableNode;
    
    /* Interpolation nodes corresponding to interpolable pts of this block */
    FldArrayI _interpolationNode;
    
    /* Indices of points requiring stringing (with no particular order) */
    FldArrayI _pointsForStrings;
    
    /* Indices of points ordered in strings */
    std::list<FldArrayI*> _strings;

    /* Vector of potentially overlapping blocks */
    std::vector<StructBlock*> _overlappingBlks;

    /* Vector of non overlapping blocks */
    std::vector<StructBlock*> _matchingBlks;

    /* Number of mesh points (optimization) */
    E_Int _nMeshPts;

    /* Bounding box of mesh */
    E_Float _xmin, _ymin, _zmin, _xmax, _ymax, _zmax;

    /* Bounding box of each cells */
    FldArrayF _bbCell;

    /* Indices of points that could be merged */
    FldArrayI _ibnd;
    
    /* Indices of the reference point of each point
       if nodes 1, 2, 4 are degenerated, ie have the same (x,y,z),
       their indice in _idg will be 1 and 1 will be the reference node
       if node 3 is not degenerated, its indice will be 3 */
    FldArrayI _idg;

    /* Unstructured triangular connectivity */
    FldArrayI _unsConnectEN;

    /* Tolerance for matching boundary condition check */
    E_Float _matchTol;
    /* Tolerance for overlap check */
    E_Float _overlapTol;
    
};
# undef FldArrayF
# undef FldArrayI
# undef FldArrayIS
# undef FldArrayB

#endif

