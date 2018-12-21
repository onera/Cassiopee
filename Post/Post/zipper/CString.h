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
#ifndef _POST_ZIPPER_CSTRING_H_
#define _POST_ZIPPER_CSTRING_H_

# include "kcore.h"
# include "StructBlock.h"
# include<vector>
# include<list>
# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI
# define FldArrayIS K_FLD::FldArrayIS

//=============================================================================
/* Class defining a 1D string of points on a mesh */
//=============================================================================
class CString
{
  public: 
    ///+ 1- Constructors / Destructor
    
    /** Constructor : build a string.
        IN : ind : array of indices of mesh points,
        IN : blk : block */
    CString(FldArrayI& ind, StructBlock* blk);
    
    /** Constructor : build a string by sub-stringing.
        IN : st : string to be sub-taken,
        IN : istart : first point to be considered,
        IN : iend : last point to be considered (included)
    */
    CString(CString* st, E_Int istart, E_Int iend);
      
    /** Destructor */
    ~CString();
    
    /** Return the block to which belongs the string */
    StructBlock* getBlock();
    
    /** Return the ind array of pts of the string */
    FldArrayI& getIndArray();
    
    /** Return the DejaVu array */
    FldArrayIS& getDejaVu();
      
    /** Get the flag array: 
     1: point interior of a segment of the string 
     0: border of a segment of the string or point not set in a segment */
    FldArrayIS& getFlag();

    /**  Return matchInfo array */
    FldArrayI& getMatchInfo();

    /** Find first non dejaVu point */
    E_Int findFirstNonDejaVu();
    
    /** Find first non flagged point */
    E_Int findFirstNonFlagged();

    /** Returns the list of potentially matching strings */
    std::vector<CString*>& getListOfMatchingStrings();
    
    /** Search for matching segments of strings and store them in the
     vector of segPairs */
    void searchForMatchingStrings( std::vector<CString*>& strings,
                                   std::vector<CString*>& strOut1,
                                   std::vector<CString*>& strOut2);

    /** Select the strings that come from overlapping blocks
        in: noStr : number of the current string in the strings vector */
    void selectStrings(std::vector<CString*>& strings, E_Int noStr);

    /** Identify segments of strings that are not set in segment pair.
        And store them in a list of arrays containing indices of points
        to be pocketted
    */
    void compRemainingSegments(std::list<FldArrayI*>& remainingSeg);

    /* Get the original string from which derives the current string */
    CString*& getOriginalString();

  private:
   
    /* Given two strings s1 and s2.
       Return the index that matches s1 on s2.
       Can return -1 : no match possible. 
       out: istart1 : first pt for segment of str1
            iend1   : last pt for segment of  str1
            istart2 : node corresponding to istart1 for str2
            iend2   : node corresponding to iend1 for str2*/
    void matchStrings(CString* s1, CString* s2,
                      std::vector<CString*>& strings,
                      E_Int& istart1, E_Int& iend1,
                      E_Int& istart2, E_Int& iend2);
    
    /* Return the nearest point of ind on another string.
       ifirst1 : in : number of pt in the _ind array
       strings : in : vector of strings
       ifirst2 : out: number of pt in the ind array nearest of ind of stringOut
       stringOut : out : string containing the nearest point */
    E_Boolean nearestPoint(E_Int ifirst1, std::vector<CString*>& strings,
                           E_Int& ifirst2, CString*& stringOut);

    /* Reset the deja Vu field to 0 */
    void resetDejaVu();
    
    /* Set dejaVu for pts of str1 between is1 and ie1 
     and for pts of str2 between is2 and ie2  */
    void setDejaVu(CString* str1, E_Int is1, E_Int ie1,
                   CString* str2, E_Int is2, E_Int ie2);

    /* Find the matching pieces of s1 and s2 */
    E_Boolean compMatching(CString* s1, CString* s2, 
                           std::vector<CString*>& strings,
                           E_Int& istart1, E_Int& istart2,
                           E_Int& iend1, E_Int& iend2);
    
    /* Test if selected segments are overlapping */
    E_Boolean areSegmentsOverlapping(E_Int is, E_Int ie,
                                     E_Int js, E_Int je);

    /* Select the most relevant segment : closest segments extremities */
    void selectClosestSegment(CString* s1, CString* s2,
                              E_Int is1, E_Int js1,
                              E_Int ie1, E_Int je1,
                              E_Int is2, E_Int js2,
                              E_Int ie2, E_Int je2,
                              E_Int& fs1, E_Int& fe1,
                              E_Int& fs2, E_Int& fe2);     

    /** For test only write outputs */
    void writeCoordOfStrings(CString* str2, E_Int is1, E_Int ie1,
                             E_Int is2, E_Int ie2);

    /* Test if extremities of strings str1 and str2 match
       and if it is the case, store information in the _matchInfo array. */
    void storeMatchingStringsInfo(std::vector<CString*>& strings,
                                  E_Int noStr1,
                                  E_Int noStr2);
    /* Find closest point of ibeg from s1 on s2
       between bounds min and max of s2
       Return the found index */
    E_Int findClosestPoint(CString* s1, CString* s2, E_Int ibeg,
                           E_Int min, E_Int max);  

  private:
    // Original string: string from which the current string comes from
    CString* _origString;
    // Indices of points forming string
    FldArrayI _ind;
    
    // Block to which the string belong
    StructBlock* _blk;
    
    // Deja vu : keeps track of already treated points of string
    // 0 : pt in not already matched with another string
    // 1 : pt already matched but on string extremity
    // 2 : pt already matched internal or 
    FldArrayIS _dejaVu;

    // Each element of the flag array :
    // 0 : point not set in any segment pair element
    // 1 : extremity point of a segment pair sample (to be assembled in 
    // a segment pair ) 
    // 2 : interior point of a segment pair sample (to be assembled in 
    // a segment pair ) 
    FldArrayIS _flag;

    // list of the strings that are candidates for matching. They come from 
    // an overlapping block of _blk
    std::vector<CString*> _listOfMatchingStrings;    

    /* Information about extremities of the current string.
       If _ind[0] or _ind[max] matches with a point of another string 
       then store in the array matchInfo
       matchInfo(0,1) : number of the string s2 (in the global strings vector)
                        that matches with pt _ind[0]
       matchInfo(0,2) : 0 if the matching pt of s2 is the minimum (ind2[0])
                        1 if the matching pt of s2 is the maximum (ind2[max2] )
       matchInfo(1,1) : number of the string s3 that matches with pt_ind[max]
       matchInfo(1,2) : 0 if the matching pt of s3 is the minimum (ind3[0])
                        1 if the matching pt of s3 is the maximum (ind3[max3] )
       par defaut (no matching) : all values are set to -1                 
    */
    FldArrayI _matchInfo;

    /* List of segments that are not already triangulated*/
    std::list<FldArrayI*> _remainingSegments;
};

# undef FldArrayF
# undef FldArrayI
# undef FldArrayIS

#endif

