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
# ifndef _POST_ZIPPER_SINGLE_SEGMENT_H_
# define _POST_ZIPPER_SINGLE_SEGMENT_H_

# include "kcore.h"
# include "StructBlock.h"
# include "SingleSegment.h"

# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI

//=============================================================================
/* Class defining a list of segments that have not been set in pairs.  */
//=============================================================================
class SingleSegment
{
  public: 
    ///+ 1- Constructors / Destructor
    
    /** Constructor */
    SingleSegment(FldArrayI& indir, StructBlock* blk);
    SingleSegment(FldArrayF& field);

    /** Destructor */
    ~SingleSegment();

    /** Return the coordinates+ cfd field of points in the segment */
    FldArrayF& getGlobalField();

  private:
    /* Coordinates+cfd field of points of the segments- already sorted */
    FldArrayF _field;
# undef FldArrayF
# undef FldArrayI
};

#endif
