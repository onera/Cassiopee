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
#ifndef _INT_TREE_NODE_H_
#define _INT_TREE_NODE_H_

namespace K_INTERP
{
// ============================================================================
/*
   The node of the ADT (binary tree of cells).*/
//=============================================================================

class IntTreeNode
{
  public:
    IntTreeNode();
    IntTreeNode(E_Int ind, E_Float xmax, E_Float ymax, E_Float zmax,
                   E_Float xmin, E_Float ymin, E_Float zmin);
    ~IntTreeNode();
    
    void setLeftAt(IntTreeNode*);
    void setRightAt(IntTreeNode*);
    IntTreeNode* getLeft() const;
    IntTreeNode* getRight() const;
    IntTreeNode** getLeftPt();
    IntTreeNode** getRightPt();
    void getCellBB(E_Float& xmax, E_Float& ymax, E_Float& zmax,
                   E_Float& xmin, E_Float& ymin, E_Float& zmin) const;
    E_Int getIndex() const;

  public: // for speed
    E_Int _ind;
    E_Float _BB[6];
    IntTreeNode* _left;
    IntTreeNode* _right;
};

// ============================================================================
// INLINE
inline void IntTreeNode::setLeftAt(IntTreeNode* left)
{
   _left = left;
}

//=============================================================================
inline void IntTreeNode::setRightAt(IntTreeNode* right)
{
   _right = right;
}

//=============================================================================
inline IntTreeNode* IntTreeNode::getLeft() const
{
   return _left;
}

//=============================================================================
inline IntTreeNode* IntTreeNode::getRight() const
{
   return _right;
}

//=============================================================================
inline IntTreeNode** IntTreeNode::getLeftPt()
{
   return &_left;
}

//=============================================================================
inline IntTreeNode** IntTreeNode::getRightPt()
{
   return &_right;
}

//=============================================================================
inline void IntTreeNode::getCellBB(E_Float& xmax, E_Float& ymax,
                                   E_Float& zmax, E_Float& xmin,
                                   E_Float& ymin, E_Float& zmin) const
{
  xmax = _BB[0];
  ymax = _BB[1];
  zmax = _BB[2];
  xmin = _BB[3];
  ymin = _BB[4];
  zmin = _BB[5];
}

//=============================================================================
inline E_Int IntTreeNode::getIndex() const
{
  return _ind;
}
// ============================================================================

}
#endif
